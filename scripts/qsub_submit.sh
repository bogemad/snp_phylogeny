#PBS -N snp_phylogeny
#PBS -l ncpus=xxxcpusxxx
#PBS -l mem=xxxmemxxx
#PBS -l walltime=200:00:00
#PBS -q workq
#PBS -o xxxbasepathxxx/logs/snp_phylogeny.out
#PBS -e xxxbasepathxxx/logs/snp_phylogeny.err

cd xxxbasepathxxx

./run.sh xxxcpusxxx &> xxxbasepathxxx/logs/snp_phylogeny.log
