#PBS -N xxxshortjobnamexxx
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=10:00:00
#PBS -q smallq
#PBS -o xxxbasepathxxx/logs/xxxshortjobnamexxx.out
#PBS -e xxxbasepathxxx/logs/xxxshortjobnamexxx.err

cd xxxbasepathxxx

scripts/multisnippy.sh xxxfilexxx xxxoutdirxxx xxxreferencexxx &> xxxbasepathxxx/logs/xxxshortjobnamexxx.log
