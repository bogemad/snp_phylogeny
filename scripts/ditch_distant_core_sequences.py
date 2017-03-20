#!/usr/bin/env python

import sys, os, subprocess, shutil

def ditch_distant_sequences(core_data):
	core_data_handle = open(core_data)
	log = open("../excluded_sequences/distant_sequences_removed_from_core_alignment.txt", 'w')
	for line in core_data_handle:
		if line.startswith("ID"):
			continue
		line_data = line.strip().split("\t")
		if float(line_data[3]) < int(sys.argv[1]):
			archive = "{}.alignment.tar.gz".format(line_data[0])
			reads_file = "../raw_data/reads/{}.fastq.gz".format(line_data[0])
			moved_reads_file = "../excluded_sequences/poor_ref_alignment/{}.fastq.gz".format(line_data[0])
			moved_archive = "../excluded_sequences/poor_ref_alignment/{}.tar.gz"".format(line_data[0])
			outline = "Sample {} coverage is too low ({}%), removing from core alignment. Reads and reference mapping data will be retained in archive: {}".format(line_data[0],line_data[3],archive)
			print(outline)
			log.write(outline + '\n')
			subprocess.check_output(["tar", "cvzf", archive, line_data[0]])
			os.rename(reads_file, moved_reads_file)
			os.rename(archive, moved_archive)
			shutil.rmtree(line_data[0])
	core_data_handle.close()
	alignment_files = [item for item in os.listdir(".") if os.path.isfile(item)]
	log.close()
	for file in alignment_files:
		if file.startswith("core"):
			os.remove(file)


ditch_distant_sequences("core.txt")

