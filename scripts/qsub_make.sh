#PBS -N make_sp
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=6:00:00
#PBS -q smallq
#PBS -o /dev/null
#PBS -e /dev/null

cd xxxbasepathxxx

make clean && make &> logs/make.log
