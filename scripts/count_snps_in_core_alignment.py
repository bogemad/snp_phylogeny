#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from Bio import AlignIO
from collections import defaultdict
import operator

# import time
import multiprocessing
# from timeit import default_timer as timer


# def time_count_subject_snps(alignment,chunksize, i, tests,threads):
	# start = timer()
	# count_subject_snps(alignment[i+tests],alignment[i],chunksize,threads)
	# end = timer()
	# return (end - start) 


# def id_peak(chunksize, chunksize_increment, alignment,threads):
	# last_average_time = float('inf')
	# tests = min(len(alignment)/2,threads)
	# while True:
		# chunksize = int(chunksize)
		# print("Attempting chunksize = %d" % chunksize)
		# time_list = [time_count_subject_snps(alignment,chunksize, i, tests,threads) for i in range(tests)]
		# average_time = sum(time_list) / float(len(time_list))
		# if average_time > last_average_time:
			# return int(chunksize-(chunksize_increment*2)), average_time
		# else:
			# chunksize += chunksize_increment
			# last_average_time = average_time


# def optimise_chunksize(alignment, threads):
	# opt_chunksize, last_average_time = id_peak(200, 200, alignment, threads)
	# opt_chunksize, last_average_time = id_peak(opt_chunksize+50, 50, alignment, threads)
	# print("Using chunksize = %d" % opt_chunksize)
	# return int(opt_chunksize+50)


# def count_chunk_snps(item):
	# subject,query = item
	# count = 0
	# if query != subject:
		# for k, q_base in enumerate(query):
			# if q_base in ["A","T","C","G"]:
				# s_base = subject[k]
				# if s_base in ["A","T","C","G"]:
					# if q_base != s_base:
						# count += 1
	# return count


# def count_subject_snps(subject,query,chunksize,threads):
	# chunk_list = [(subject[i:i+chunksize].seq,query[i:i+chunksize].seq) for i in range(0, len(query), chunksize)]
	# count_list = []
	# p = multiprocessing.Pool(threads)
	# count_list = p.map_async(count_chunk_snps,chunk_list,chunksize=1).get()
	# p.close()
	# p.join()
	# return sum(count_list)


# def main():
	# alignment_file = sys.argv[1]
	# threads = int(sys.argv[3])
	# alignment = AlignIO.read(alignment_file,"fasta")
	# print("Optimizing script...")
	# chunksize = optimise_chunksize(alignment, threads)
	# column_names = [x.id for x in alignment]
	# df = pd.DataFrame(data=np.zeros((0,len(column_names))), columns=column_names)
	# for query in alignment:
		# snp_count_d = {}
		# for subject in alignment:
			# snp_count_d[subject.id] = count_subject_snps(subject,query,chunksize,threads)
		# print("Adding %s results to dataframe" % query.id)
		# df_1 = pd.DataFrame(snp_count_d, index=[query.id])
		# df = df.append(df_1)
	# df.to_csv(sys.argv[2], sep="\t")

def gen_chunks(sequence, threads):
	chunk_num = int(threads)
	chunk_size = int(len(sequence[0])/chunk_num)+1
	return [(i, sequence[:,i:i + chunk_size]) for i in range(0, len(sequence[0]), chunk_size)]


def import_array_in_chunks(indexed_chunk):
	index, aln = indexed_chunk
	chunk = np.array([list(rec) for rec in aln], np.character)
	indices = [i for i, v in enumerate(chunk[0]) if (b'N' in chunk[:,i]) or (b'-' in chunk[:,i])]
	chunk = np.delete(chunk, indices, axis=1)
	return index, chunk

def multi_import(chunk_list, threads):
	p = multiprocessing.Pool(threads)
	results = p.map_async(import_array_in_chunks,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	return results


def import_array(aln, threads):
	chunk_list = gen_chunks(aln, threads)
	array_chunk_list = multi_import(chunk_list, threads)
	return reassemble_array(array_chunk_list)


def reassemble_array(array_chunk_list):
	l = [array_chunk for i, array_chunk in sorted(array_chunk_list,key=operator.itemgetter(0))]
	return np.concatenate(l , axis=1)


def multijob(input):
	query_id, subject_id, query, subject = input
	snps = int(np.count_nonzero(query != subject))
	print("{} SNPs found between isolates {} and {}".format(snps, query_id, subject_id))
	return (query_id, subject_id, snps)


def multirun(input_list, threads):
	p = multiprocessing.Pool(threads)
	result = p.map_async(multijob,input_list,chunksize=1)
	prev_jobs_remaining = len(input_list)
	while not result.ready():
		jobs_remaining = result._number_left
		if jobs_remaining != prev_jobs_remaining:
			print("{} of {} sequence pairs remaining to be compaired".format(result._number_left, len(input_list)))
		prev_jobs_remaining = jobs_remaining
	results_list = result.get(999999999999999)
	p.close()
	p.join()
	return results_list


def gen_input_list(array, id_list):
	input_list = []
	for i, query in enumerate(array):
		for j, subject in enumerate(array):
			input_list.append((id_list[i], id_list[j], query, subject))
	return input_list

def sort_results(results):
	snp_count_d = defaultdict(dict)
	for query_id, subject_id, sum in results:
		snp_count_d[query_id][subject_id] = sum
	return snp_count_d


def generate_data_frame(snp_count_d,id_list):
	df = pd.DataFrame(data=np.zeros((0,len(id_list))), columns=id_list)
	for query_id in sorted(id_list):
		print("Adding %s results to dataframe" % query_id)
		df_1 = pd.DataFrame(snp_count_d[query_id], index=[query_id])
		df = df.append(df_1)
	return df


def validate_multi_import(array, array2):
	for i in range(len(array)):
		if multijob((i ,i, array[i], array2[i]))[2] != 0:
			print("Check multiprocessing array input is working")
			sys.exit(1)
	print("Multiprocessing array input is working fine")


def main():
	alignment_file = sys.argv[1]
	threads = int(sys.argv[3])
	print("Importing alignment...")
	aln = AlignIO.read(alignment_file,"fasta")
	print("Cleaning alignment of gaps...")
	array = import_array(aln, threads)
	# array2 = np.array([list(rec) for rec in aln], np.character)
	# validate_multi_import(array, array2)
	id_list = [x.id for x in aln]
	input_list = gen_input_list(array, id_list)
	results = multirun(input_list, threads)
	print("Sorting results...")
	snp_count_d = sort_results(results)
	df = generate_data_frame(snp_count_d,id_list)
	df.to_csv(sys.argv[2], sep="\t")
	print("Done.")

if __name__ == '__main__':
	main()
