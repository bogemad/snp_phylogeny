#!/usr/bin/env python

import sys, os, shutil, subprocess, multiprocessing


def gen_reads_dict(file):
	with open(file) as file_handle:
		reads_dict = {}
		for line in file_handle:
			if line.startswith('Uberstrain'):
				continue
			data = line.strip().split('\t')
			name = data[1].replace(" ","_").replace("(","").replace(")","")
			reads_dict[name] = [sub_data.split(";")[0].strip('"') for sub_data in data[2].split(",")]
	return reads_dict


def download_reads(reads_list, name):
	base_path = os.path.dirname(os.path.dirname(__file__))
	temp = os.path.join(base_path,".temp")
	if os.path.isdir(temp) == False:
		os.mkdir(temp)
	print("%s: Downloading reads..." % name)
	for SRR in reads_list:
		if os.path.exists((SRR + ".fastq.gz")) == False:
			subprocess.call(["fastq-dump", "--gzip", "--skip-technical", "--read-filter", "pass", "--dumpbase", "--split-spot", "--clip", "--outdir", temp, SRR])
	print("%s: Read download complete." % name)
	return [os.path.join(temp,(SRR + "_pass.fastq.gz")) for SRR in reads_list]


def merge_reads(name, filenames):
	base_path = os.path.dirname(os.path.dirname(__file__))
	data_dir = os.path.join(base_path,"raw_data","reads")
	merged_fastq = os.path.join(data_dir,(name + ".fastq.gz"))
	print("%s: Merging reads..." % name)
	if os.path.exists(merged_fastq) == False:
		with open(merged_fastq, 'wb') as outfile:
			for fname in filenames:
				with open(fname, 'rb') as infile:
					shutil.copyfileobj(infile, outfile)
				os.remove(fname)
	print("%s: Read merge complete." % name)
	return merged_fastq


def download_align_reads(data):
	name, reads_list = data
	filenames = download_reads(reads_list, name)
	merged_fastq = merge_reads(name, filenames)


if __name__ == "__main__":
	print("Generating download list...")
	SRR_dict = gen_reads_dict(sys.argv[1])
	threads = sys.argv[2]
	# for item in SRR_dict.items():
		# download_align_reads(item)
	p = multiprocessing.Pool(int(threads))
	p.map_async(download_align_reads, SRR_dict.items(), chunksize=1).get()
	p.close()
	p.join()
	print("Done.")

