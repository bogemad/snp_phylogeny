#!/bin/bash

file=$1
outdir=$2
reference=$3
conda_bin_path=`readlink -f .mc/bin`
scripts_path=`readlink -f scripts`
PATH=$conda_bin_path:$scripts_path:$PATH

name=`get_samplenames.py $file`

if [ $name == 'fail' ]
then
	echo "Failed to determine samplename for file path: $file. Check if this is a fastq or fastq.gz file."
	exit 1
fi

read_format=`id_read_format.py $file`

if [ $read_format == "fastq" ]; do
	snippy --prefix $name --cpus 1 --outdir $outdir/$name --ref $reference --peil $file || (echo "Alignment failed! Check if $file is correct fastq file." && exit 1)
elif [ $read_format == "fasta" ]; do
	snippy --prefix $name --cpus 1 --outdir $outdir/$name --ref $reference --ctgs $file || (echo "Alignment failed! Check if $file is correct fasta file." && exit 1)
else
	echo "Cannot determine format of read file: $file. 
	Please ensure it is either fastq (for raw sequencing reads) or fasta (for assembled genomes)"
fi
