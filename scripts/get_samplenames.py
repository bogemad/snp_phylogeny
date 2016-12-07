#!/usr/bin/env python

import os, sys

filename = sys.argv[1]
basename = os.path.basename(filename)

name, ext = os.path.splitext(basename)

if ext == ".gz":
	name, ext = os.path.splitext(name)


if ext == ".fastq" or ext == ".fq":
	name2, ext2 = os.path.splitext(name)
	if ext2 == ".lf":
		name = name2
	print(name)
else:
	print("fail")