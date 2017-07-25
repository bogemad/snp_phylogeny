#!/bin/bash

base_path=$PWD
echo $base_path

cp -av $base_path/test/raw_data/reads/* $base_path/raw_data/reads/
cp -av $base_path/test/raw_data/reference_sequence/* $base_path/raw_data/reference_sequence/

./run.sh 4 90

rm -rf $base_path/raw_data/reads/* $base_path/raw_data/reference_sequence/* $base_path/analysis_results/* $base_path/excluded_sequences
