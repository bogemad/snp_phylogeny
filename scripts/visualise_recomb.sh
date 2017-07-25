#!/bin/bash

source activate gubbins-env

python ../scripts/gubbins_drawer.py -o filtered_core_aln.recombination_predictions.pdf -t filtered_core_aln.final_tree.tre filtered_core_aln.recombination_predictions.embl

source deactivate gubbins-env
