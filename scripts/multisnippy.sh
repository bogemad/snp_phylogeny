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

snippy --prefix $name --cpus 1 --outdir $outdir/$name --ref $reference --peil $file || (echo "Alignment failed! Check if $filename is correct fastq file." && exit 1)
