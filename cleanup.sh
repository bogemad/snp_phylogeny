#!/bin/bash

base_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
temp=$base_path/.temp
logs=$base_path/logs

rm -rf $temp
mkdir -p $temp
mkdir -p $logs

cd $base_path

mv excluded_sequences/poor_ref_alignment/*.tar.gz analysis_results
mv excluded_sequences/poor_ref_alignment/* raw_data/reads

cd analysis_results

parallel 'tar xvzf {}' ::: `ls *.tar.gz`

rm -rf *.tar.gz core* RAxML* total.snp*
cd ..
rm -rf excluded_sequences

for dir in `ls analysis_results`; do
if `ls raw_data/reads/${dir}* &> /dev/null`; then
echo "$dir has a read file and will be kept"
else
echo "$dir has no read file and will to be deleted"
rm -rf analysis_results/$dir
fi
done
