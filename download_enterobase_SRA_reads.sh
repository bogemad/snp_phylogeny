#!/bin/bash

conda_bin_path=`readlink -f .mc/bin`
scripts_path=`readlink -f scripts`

PATH=$conda_bin_path:$scripts_path:$PATH

download_enterobase_SRA_reads.py $1 $2
