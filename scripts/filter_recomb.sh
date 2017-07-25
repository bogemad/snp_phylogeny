#!/bin/bash

threads=$1
PATH=$conda_bin_path:$scripts_path:$PATH

echo "Running recombination filter..."

source activate gubbins-env
run_gubbins.py -v -i 10 -p filtered_core_aln -c $threads core.full.trimmed.aln || gubbins_fail=true

if [ "$gubbins_fail" == true ]; then
	echo "gubbins using RAxML only method has failed. Retrying with fastree for first iteration."
	gubbins_fail=false
	rm -rf core.full.trimmed.aln.*
	run_gubbins.py -v -i 10 --tree_builder hybrid -p filtered_core_aln -c $threads core.full.trimmed.aln || gubbins_fail=true
fi

if [ "$gubbins_fail" == true ]; then
	echo "gubbins using hybrid RAxML/FastTree method has failed. Retrying with fastree for all iterations."
	gubbins_fail=false
	rm -rf core.full.trimmed.aln.*
	run_gubbins.py -v -i 10 --tree_builder fasttree -p filtered_core_aln -c $threads core.full.trimmed.aln || (echo "Running gubbins using all methods have failed. Please examine your alignment or consider removing highly divergent sequences. Looking at total.snp_counts.stats is a good place to start. Additionally consider using a different reference sequence."; exit 1)
fi

source deactivate gubbins-env
