#!/bin/bash
#PBS -N snp_phy_cleanup
#PBS -l ncpus=4
#PBS -l mem=4gb
#PBS -l walltime=2:00:00
#PBS -q workq
#PBS -o xxxbasepathxxx/logs/snp_phy_cleanup.out
#PBS -e xxxbasepathxxx/logs/snp_phy_cleanup.err

cd xxxbasepathxxx

mv excluded_sequences/poor_ref_alignment/*.tar.gz analysis_results
mv excluded_sequences/poor_ref_alignment/* raw_data/reads

cd analysis_results

parallel -j 4 'tar xvzf {}' ::: `ls *.tar.gz`

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
