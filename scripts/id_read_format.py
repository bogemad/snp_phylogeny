#!/usr/bin/env python

import sys, os, gzip

def main():
	file = sys.argv[1]
	if os.path.splitext(file)[1] == ".gz":
		handle = gzip.open(file,'rt')
	else:
		handle = open(file)
	for line in handle:
		if line.strip() == "":
			continue
		if line.strip().startswith('@'):
			print("fastq")
		elif line.strip().startswith(">"):
			print("fasta")
		else:
			sys.exit(1)
		return

main()