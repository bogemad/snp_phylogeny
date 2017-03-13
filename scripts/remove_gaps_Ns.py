#!/usr/bin/env python

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import multiprocessing
import sys
import numpy as np
import operator


def gen_chunks(sequence, threads):
	chunk_num = int(threads)
	chunk_size = int(len(sequence[0])/chunk_num)+1
	return [(i, sequence[:,i:i + chunk_size]) for i in range(0, len(sequence[0]), chunk_size)]


def import_edit_array_in_chunks(indexed_chunk):
	index, aln = indexed_chunk
	chunk = np.array([list(rec) for rec in aln], np.character)
	indices = [i for i, v in enumerate(chunk[0]) if (b'N' in chunk[:,i]) or (b'-' in chunk[:,i])]
	chunk = np.delete(chunk, indices, axis=1)
	return index, chunk


def multi_import(chunk_list, threads):
	p = multiprocessing.Pool(threads)
	results = p.map_async(import_edit_array_in_chunks,chunk_list,chunksize=1).get()
	p.close()
	p.join()
	return results


def filter_gaps_Ns(aln, threads):
	chunk_list = gen_chunks(aln, threads)
	array_chunk_list = multi_import(chunk_list, threads)
	return reassemble_array(array_chunk_list)


def reassemble_array(array_chunk_list):
	l = [array_chunk for i, array_chunk in sorted(array_chunk_list,key=operator.itemgetter(0))]
	return np.concatenate(l , axis=1)


def reconstruct_alignment(array,aln):
	for i, rec in enumerate(aln):
		sequence_bytes = b''.join(array[i])
		rec.seq = Seq(str(sequence_bytes,'utf-8'),IUPAC.unambiguous_dna)
	return aln



def main():
	print("Importing alignment...")
	aln = AlignIO.read(sys.argv[1],"fasta")
	print("Cleaning alignment...")
	array = filter_gaps_Ns(aln, int(sys.argv[2]))
	print("Finished cleaning. Exporting to file...")
	AlignIO.write(reconstruct_alignment(array,aln), sys.argv[3], "fasta")
	print("Done.")


if __name__ == "__main__":
	main()
