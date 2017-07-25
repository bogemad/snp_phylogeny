#!/bin/bash

git clone https://github.com/sanger-pathogens/gubbins.git
home_dir=`dirname $1`
source activate gubbins-env 
cd gubbins 
autoreconf -i 
./configure --prefix=${1}/envs/gubbins-env
make 
make install 
cd python 
python setup.py install
mv $home_dir/gubbins/python/scripts/gubbins_drawer.py $home_dir/scripts
cd ../..
rm -rf gubbins
source deactivate