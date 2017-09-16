#!/usr/bin/env python

import sys
import os
import gzip
import subprocess
import shutil
import re
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import multiprocessing
import numpy as np
import pandas as pd
import operator
from collections import defaultdict

def parse_sequence_file(sequence_file):
	fasta_parse = [x for x in SeqIO.parse(sequence_file,'fasta')]
	gb_parse = [x for x in SeqIO.parse(sequence_file,'gb')]
	if len(fasta_parse) > 0:
		return fasta_parse, False
	elif len(gb_parse) > 0:
		return gb_parse, True
	else:
		print("Unable to parse reference sequence file, please check if it is fasta or gb format.")
		sys.exit(1)


def write_new_sequence_file(cleaned_sequences, annotations_present, outfile):
	print("Writing sequence to %s" % outfile)
	if annotations_present == True:
		outfile = outfile + ".gbk"
		null = SeqIO.write(cleaned_sequences, outfile, "gb")
	elif annotations_present == False:
		outfile = outfile + ".fasta"
		null = SeqIO.write(cleaned_sequences, outfile, "fasta")


def check_degenerate(base):
	if base.upper() in ['A','T','G','C','N','-']:
		return False
	return True

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[:, i:i + n]

def gen_ref_chunks(sequence, threads):
	chunk_num = threads ** 2
	chunk_size = int(len(sequence)/chunk_num)
	if chunk_size < 1:
		chunk_size = 1
	return [(i, sequence[i:i + chunk_size]) for i in range(0, len(sequence), chunk_size)]

def filter_degenerates(indexed_sequence_chunk):
	index, sequence = indexed_sequence_chunk
	degenerate_locations = [index+i for i, base in enumerate(sequence) if check_degenerate(base) == True]
	return index, degenerate_locations

def clean_sequence(results_list, sequence):
	result_dict = {i:j for i,j in results_list}
	keys = list(result_dict.keys())
	degenerate_locations = []
	for i in sorted(keys):
		degenerate_locations += result_dict[i]
	print("Total %d degenerate bases removed from sequence" % len(degenerate_locations))
	mutable_seq = sequence.seq.tomutable()
	for location in degenerate_locations:
		mutable_seq[location] = "N"
	new_seq = mutable_seq.toseq()
	sequence.seq = new_seq
	return sequence

def id_read_format(file):
	if os.path.splitext(file)[1] == ".gz":
		handle = gzip.open(file,'rt')
	else:
		handle = open(file)
	for line in handle:
		if line.strip() == "":
			continue
		if line.strip().startswith('@'):
			return "fastq"
		elif line.strip().startswith(">"):
			return "fasta"
		else:
			print("Can't properly determine read format for file {}. Please check your read files are in fastq for raw sequencing reads) or fasta (for assembled genomes) format".format(os.path.basename(file)))
			sys.exit(1)


def get_samplename(filename):
	basename = os.path.basename(filename)
	name, ext = os.path.splitext(basename)
	if ext == ".gz":
		name, ext = os.path.splitext(name)
	if ext == ".fastq" or ext == ".fq" or ext == ".fasta" or ext == ".fa" or ext == ".fna" or ext == ".fsa":
		name2, ext2 = os.path.splitext(name)
		if ext2 == ".lf":
			name = name2
		return name
	else:
		print("Failed to determine samplename for file: {}. Check if this is a fastq or fastq.gz file.".format(basename))
		sys.exit(1)


def submit_hpc_script(temp_dir, read_path, outdir, reference, threads, mem, queue, base_path, id_threshold):
	with open(os.path.join(temp_dir,"qsub.sh"), 'w') as sub_file:
		sub_file.write("#PBS -N {}\n".format(os.path.basename(read_path)[:15]))
		sub_file.write("#PBS -l ncpus={}\n".format(threads))
		sub_file.write("#PBS -l mem={}\n".format(mem))
		sub_file.write("#PBS -l walltime=2:00:00\n")
		sub_file.write("#PBS -q {}\n".format(queue))
		sub_file.write("#PBS -o {}/logs/snp_phylogeny.out\n".format(base_path))
		sub_file.write("#PBS -e {}/logs/snp_phylogeny.err\n\n".format(base_path))
		sub_file.write("export PATH={}/.mc/bin:$PATH\n".format(base_path))
		sub_file.write("cd {}\n".format(base_path))
		sub_file.write("python {0}/scripts/run.py {0} {1} {2} {3} {4} true {5} {6} {8} &> {0}/logs/{7}_snippy.log\n".format(base_path, os.path.dirname(read_path), reference, outdir, threads, id_threshold, read_path, get_samplename(read_path), mem))
	p = subprocess.Popen(["qsub", "-V", os.path.join(temp_dir,"qsub.sh")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output = p.communicate()[0].decode('utf-8')
	id = output.strip()
	os.remove(os.path.join(temp_dir,"qsub.sh"))
	return id


def submit_tree_hpc_script(temp_dir, read_path, outdir, reference, threads, mem, queue, base_path, id_threshold, alignment_ids):
	with open(os.path.join(temp_dir,"qsub.sh"), 'w') as sub_file:
		sub_file.write("#PBS -N sp_tree_build\n")
		sub_file.write("#PBS -l ncpus={}\n".format(threads))
		sub_file.write("#PBS -l mem={}\n".format(mem))
		sub_file.write("#PBS -l walltime=200:00:00\n")
		sub_file.write("#PBS -q {}\n".format(queue))
		sub_file.write("#PBS -o {}/logs/snp_phylogeny.out\n".format(base_path))
		sub_file.write("#PBS -e {}/logs/snp_phylogeny.err\n\n".format(base_path))
		sub_file.write("export PATH={}/.mc/bin:$PATH\n".format(base_path))
		sub_file.write("cd {}\n".format(base_path))
		sub_file.write("python {0}/scripts/run.py {0} {1} {2} {3} {4} true {5} {6} {7} &>> {0}/logs/snp_phylogeny.log\n".format(base_path, os.path.dirname(read_path), reference, outdir, threads, id_threshold, 'all', mem))
	subprocess.run(["qsub", "-V", "-W", "depend=afterok:{}".format(":".join(alignment_ids)), os.path.join(temp_dir,"qsub.sh")], check=True, stderr=subprocess.STDOUT)
	os.remove(os.path.join(temp_dir,"qsub.sh"))



def ditch_distant_sequences(core_data, id_threshold):
	core_data_handle = open(core_data)
	log = open("../excluded_sequences/distant_sequences_removed_from_core_alignment.txt", 'w')
	bad_seqs = []
	for line in core_data_handle:
		if line.startswith("ID"):
			continue
		line_data = line.strip().split("\t")
		if float(line_data[3]) < int(id_threshold):
			archive = "{}.alignment.tar.gz".format(line_data[0])
			reads_file = "../raw_data/reads/{}".format(find_source_file(line_data[0]))
			moved_reads_file = "../excluded_sequences/poor_ref_alignment/{}".format(find_source_file(line_data[0]))
			moved_archive = "../excluded_sequences/poor_ref_alignment/{}.tar.gz".format(line_data[0])
			outline = "Sample {} coverage is too low ({}%), removing from core alignment. Reads and reference mapping data will be retained in archive: {}".format(line_data[0],line_data[3],archive)
			print(outline)
			log.write(outline + '\n')
			subprocess.check_output(["tar", "cvzf", archive, line_data[0]])
			os.rename(reads_file, moved_reads_file)
			os.rename(archive, moved_archive)
			shutil.rmtree(line_data[0])
			bad_seqs.append(line_data[0])
	core_data_handle.close()
	alignment_files = [item for item in os.listdir(".") if os.path.isfile(item)]
	log.close()
	for file in alignment_files:
		if file.startswith("core"):
			os.remove(file)
	return bad_seqs


def find_source_file(name):
	reads_files = os.listdir("../raw_data/reads")
	for file in reads_files:
		searchObj = re.search(name,file)
		if searchObj:
			return file



def gen_chunks(sequence, threads):
	chunk_num = int(threads)
	chunk_size = int(len(sequence[0])/chunk_num)+1
	return [(i, sequence[:,i:i + chunk_size]) for i in range(0, len(sequence[0]), chunk_size)]


def import_edit_array_in_chunks(indexed_chunk):
	index, aln = indexed_chunk
	chunk = np.array([list(rec) for rec in aln], np.character)
	chunk[chunk == b'N'] = b'-'
	return index, chunk

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


def multi_edit_import(chunk_list, threads):
	p = multiprocessing.Pool(threads)
	results = p.map_async(import_edit_array_in_chunks,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	return results

def import_array(aln, threads):
	chunk_list = gen_chunks(aln, threads)
	array_chunk_list = multi_import(chunk_list, threads)
	return reassemble_array(array_chunk_list)


def filter_gaps_Ns(aln, threads):
	chunk_list = gen_chunks(aln, threads)
	array_chunk_list = multi_edit_import(chunk_list, threads)
	return reassemble_array(array_chunk_list)


def reassemble_array(array_chunk_list):
	l = [array_chunk for i, array_chunk in sorted(array_chunk_list,key=operator.itemgetter(0))]
	return np.concatenate(l , axis=1)


def reconstruct_alignment(array,aln):
	for i, rec in enumerate(aln):
		sequence_bytes = b''.join(array[i])
		rec.seq = Seq(str(sequence_bytes,'utf-8'),IUPAC.unambiguous_dna)
	return aln


def replace_Ns_with_gaps(alignment_file, threads, outfile):
	print("Importing alignment...")
	aln = AlignIO.read(alignment_file,"fasta")
	print("Filtering Ns...")
	array = filter_gaps_Ns(aln, int(threads))
	print("Finished cleaning. Exporting to file...")
	AlignIO.write(reconstruct_alignment(array,aln), outfile, "fasta")
	print("Done.")

def count_snps(alignment_file, outfile, statsfile, threads):
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
	stats = df.describe()
	stats.to_csv(statsfile, sep="\t")
	df.to_csv(outfile, sep="\t")
	print("Done.")


def remove_degenerate_bases(sequence_file, outfile, threads):
	print("Scanning for degenerate reference bases...")
	sequences, annotations_present = parse_sequence_file(sequence_file)
	cleaned_sequences = []
	for sequence in sequences:
		print(len(sequence))
		chunk_list = gen_ref_chunks(sequence, threads)
		p = multiprocessing.Pool(threads)
		result = p.map_async(filter_degenerates, chunk_list, chunksize=1)
		results_list = result.get()
		p.close()
		p.join()
		cleaned_sequences.append(clean_sequence(results_list, sequence))
	write_new_sequence_file(cleaned_sequences, annotations_present, outfile)
