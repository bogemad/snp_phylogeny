#!/usr/bin/env python

import sys
from BCBio import GFF
from Bio import AlignIO

infile = sys.argv[1]

def change_fse_item(new, old, list):
	new_list = []
	for item in list:
		if item == old:
			new_list.append(new)
		else:
			new_list.append(item)
	return new_list

def gen_se_list(infile):
	with open(infile) as inhandle:
		for entry in GFF.parse(inhandle):
			list = [(int(feature.location.start),int(feature.location.end)) for feature in entry.features]
	return list

def flatten_list(infile, list):
	new_list = []
	for start, end in list:
		for fse_start, fse_end in new_list:
			if start > fse_start and end < fse_end:
				break
			if start < fse_start and end > fse_end:
				new_list = change_fse_item((start,end), (fse_start,fse_end), new_list)
				break
			if start > fse_start and start < fse_end:
				if end > fse_end:
					new_list = change_fse_item((fse_start,end), (fse_start,fse_end), new_list)
				break
			if end > fse_start and end < fse_end:
				if start < fse_start:
					new_list = change_fse_item((start,fse_end), (fse_start,fse_end), new_list)
				break
		else:
			new_list.append((start,end))
	if len(list) != len(new_list):
			return new_list, True
	return new_list, False

se_list = gen_se_list(infile)
continu = True
while continu == True:
	se_list, continu = flatten_list(infile, se_list)

alignment_file = sys.argv[2]
alignment = AlignIO.read(alignment_file, 'fasta')

key = len(alignment[0])
for i, (start, end) in enumerate(sorted(se_list)):
	if i == 0:
		print("Adding segment 0 to %d" % start)
		new_alignment = alignment[:,:start]
	elif i == len(se_list)-1:
		print("Adding segment %d to %d" % (end,len(alignment[0])))
		new_alignment += alignment[:,end:]
	else:
		next_start, null = sorted(se_list)[i+1]
		print("Adding segment %d to %d" % (end,next_start))
		new_alignment += alignment[:,end:next_start]

outfile = sys.argv[3]

count = AlignIO.write(new_alignment,outfile,"fasta")
