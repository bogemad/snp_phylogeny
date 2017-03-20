#!/usr/bin/env python

import os, sys, re, multiprocessing
from Bio import SeqIO

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

def gen_chunks(sequence, threads):
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


sequence_file = sys.argv[1]
threads = int(sys.argv[3])
print("Scanning for degenerate reference bases...")
sequences, annotations_present = parse_sequence_file(sequence_file)
cleaned_sequences = []
for sequence in sequences:
	print(len(sequence))
	chunk_list = gen_chunks(sequence, threads)
	p = multiprocessing.Pool(threads)
	result = p.map_async(filter_degenerates, chunk_list, chunksize=1)
	results_list = result.get()
	p.close()
	p.join()
	cleaned_sequences.append(clean_sequence(results_list, sequence))

write_new_sequence_file(cleaned_sequences, annotations_present, sys.argv[2])
