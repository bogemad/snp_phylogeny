#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from Bio import AlignIO
import time
import multiprocessing
from timeit import default_timer as timer


def time_count_subject_snps(alignment,chunksize, i, tests,threads):
	start = timer()
	count_subject_snps(alignment[i+tests],alignment[i],chunksize,threads)
	end = timer()
	return (end - start) 


def id_peak(chunksize, chunksize_increment, alignment,threads):
	last_average_time = float('inf')
	tests = min(len(alignment)/2,threads)
	while True:
		chunksize = int(chunksize)
		print("Attempting chunksize = %d" % chunksize)
		time_list = [time_count_subject_snps(alignment,chunksize, i, tests,threads) for i in range(tests)]
		average_time = sum(time_list) / float(len(time_list))
		if average_time > last_average_time:
			return int(chunksize-(chunksize_increment*2)), average_time
		else:
			chunksize += chunksize_increment
			last_average_time = average_time


def optimise_chunksize(alignment, threads):
	opt_chunksize, last_average_time = id_peak(200, 200, alignment, threads)
	opt_chunksize, last_average_time = id_peak(opt_chunksize+50, 50, alignment, threads)
	print("Using chunksize = %d" % opt_chunksize)
	return int(opt_chunksize+50)


def count_chunk_snps(item):
	subject,query,alignment = item
	count = 0
	if query != subject:
		for k, q_base in enumerate(query):
			if q_base in ["A","T","C","G"]:
				s_base = subject[k]
				if s_base in ["A","T","C","G"]:
					if q_base != s_base and 'N' in alignment[:,k] == False:
						count += 1
	return count


def count_subject_snps(alignment,subject,query,chunksize,threads):
	chunk_list = [(subject[i:i+chunksize].seq,query[i:i+chunksize].seq,alignment) for i in range(0, len(query), chunksize)]
	count_list = []
	p = multiprocessing.Pool(threads)
	count_list = p.map_async(count_chunk_snps,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	return sum(count_list)


def main():
	alignment_file = sys.argv[1]
	threads = int(sys.argv[3])
	alignment = [x for x in AlignIO.read(alignment_file,"fasta")]
	print("Optimizing script...")
	chunksize = optimise_chunksize(alignment, threads)
	column_names = [x.id for x in alignment]
	df = pd.DataFrame(data=np.zeros((0,len(column_names))), columns=column_names)
	for query in alignment:
		snp_count_d = {}
		for subject in alignment:
			snp_count_d[subject.id] = count_subject_snps(alignment,subject,query,chunksize,threads)
		print("Adding %s results to dataframe" % query.id)
		df_1 = pd.DataFrame(snp_count_d, index=[query.id])
		df = df.append(df_1)
	df.to_csv(sys.argv[2], sep="\t")

if __name__ == '__main__':
	main()

