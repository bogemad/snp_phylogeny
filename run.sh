#!/bin/bash

base_path=$PWD
reads_dir=`readlink -f raw_data/reads`

reference=`readlink -f raw_data/reference_sequence`
outdir=`readlink -f analysis_results`
threads=$1
conda_bin_path=`readlink -f .mc/bin`
scripts_path=`readlink -f scripts`
count_snps_in_core_alignment_fail=false
gubbins_drawer_fail=false
PATH=$conda_bin_path:$scripts_path:$PATH
temp=`readlink -f .temp`
mkdir -p $temp

if [ -z "$2" ]
then
  hpc=false
else
  hpc=$3
  mem=$4
fi
id_threshold=$2

echo ""
echo "Reads directory = $reads_dir"
echo "Reference sequence file = $reference"
echo "Output directory = $outdir"
echo "Process threads = $threads"
echo "Minimum reference coverage threshold = $id_threshold%"
echo ""

python $base_path/scripts/run.py $base_path $reads_dir $reference $outdir $threads $hpc $id_threshold all $mem

# for dir in `find $base_path/analysis_results -maxdepth 1 -mindepth 1 -type d -exec basename {} \;`; do
# if `ls $base_path/raw_data/reads/${dir}* &> /dev/null`; then
# useless_var=true
# else
# echo "$dir has no associated read file and will be deleted"
# rm -rf $base_path/analysis_results/$dir
# fi
# done


# ref_files=`find raw_data/reference_sequence -type f | wc -l`
# if [[ $ref_files > 1 ]]; then
# echo "Make sure there is only one reference file in raw_data/reference_sequence directory."
# exit 1
# fi

# if [[ "$reference" == *.gbk || "$reference" == *.gb || "$reference" == *.gbff ]]; then
  # refsuff=".gbk"
# elif [[ "$reference" == *.fasta || "$reference" == *.fa || "$reference" == *.fna ]]; then
  # refsuff=".fasta"
# else
  # echo "Can't id reference sequence format. Ensure you have a Genbank (.gb, .gbk, .gbff) or Fasta (.fasta,.fa,.fna) formatted sequence."
  # exit 1
# fi

# if [ ! -f $outdir/reference$refsuff ]; then
  # remove_degenerate_bases.py $reference $outdir/reference $threads || exit 1
# else
  # echo "Reference sequence has already been processed. Skipping this step..."
# fi

# reference=$outdir/reference$refsuff
# for file in `find $reads_dir -type f`; do
  # name=`get_samplenames.py $file`
  # snippy_done=()
  # snippy_not_done=()
  # if `ls $base_path/analysis_results/$name &> /dev/null`; then
    # snippy_done+=("$file")
  # else
    # snippy_not_done+=("$file")
  # fi
# done

# echo "${snippy_done[@]}"
# echo "${snippy_not_done[@]}"

# exit

# for reads in "${snippy_not_done[@]}"; do
	# bash multisnippy.sh $reads $outdir $reference $threads
# done

# parallel -j $threads --ungroup multisnippy.sh {} $outdir $reference ::: "${snippy_not_done[@]}"


# if [ "$exit_on_end" == "true" ]
# then
# echo "Snippy jobs have been submitted to the HPC smallq for processing. Please restart the pipeline when these jobs have been completed."
# exit
# fi

# cd $outdir
# rm $reference

# gubbins_fail=false

# if [ ! -f $outdir/core.full.aln ]; then
  # snippy-core * || exit 1
  # mkdir -p ../excluded_sequences/poor_ref_alignment
  # ditch_distant_core_sequences.py $id_threshold
  # snippy-core * || exit 1
# else
  # echo "Core alignment has already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/core.full.trimmed.aln ]; then
  # replace_Ns_with_gaps.py core.full.aln $threads core.full.trimmed.aln || exit 1
# else
  # echo "Core alignment has already been trimmed. Skipping this step..."
# fi

# if [ ! -f $outdir/total.snp_counts.txt ]; then
  # count_snps_in_core_alignment.py core.full.trimmed.aln total.snp_counts.txt total.snp_counts.stats $threads
# else
  # echo "Total snp count matrix has already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/filtered_core_aln.final_tree.tre ]; then
  # echo "Running recombination filter..."
  # source activate gubbins-env
  # run_gubbins.py -v -i 10 -p filtered_core_aln -c $threads core.full.trimmed.aln || gubbins_fail=true
  # if [ "$gubbins_fail" == true ]; then
    # echo "gubbins using RAxML only method has failed. Retrying with fastree for first iteration."
    # gubbins_fail=false
    # rm -rf core.full.trimmed.aln.*
    # run_gubbins.py -v -i 10 --tree_builder hybrid -p filtered_core_aln -c $threads core.full.trimmed.aln || gubbins_fail=true
  # fi
  # if [ "$gubbins_fail" == true ]; then
    # echo "gubbins using hybrid RAxML/FastTree method has failed. Retrying with fastree for all iterations."
    # gubbins_fail=false
    # rm -rf core.full.trimmed.aln.*
    # run_gubbins.py -v -i 10 --tree_builder fasttree -p filtered_core_aln -c $threads core.full.trimmed.aln || (echo "Running gubbins using all methods have failed. Please examine your alignment or consider removing highly divergent sequences. Looking at total.snp_counts.stats is a good place to start..."; exit 1)
  # fi
  # source deactivate gubbins-env
# else
  # echo "Recombination filtering by gubbins has already been done. Skipping this step..."
# fi

# if [ ! -f $outdir/filtered_core_aln.recombination_predictions.pdf ]; then
  # gubbins_drawer.py -o filtered_core_aln.recombination_predictions.pdf -t filtered_core_aln.final_tree.tre filtered_core_aln.recombination_predictions.embl || gubbins_drawer_fail=true
# else
  # echo "Recombination prediction visualisation has already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/filtered_core_aln.gubbins_filtered.aln ]; then
  # filter_recombinations.py filtered_core_aln.recombination_predictions.gff core.full.trimmed.aln filtered_core_aln.gubbins_filtered.aln || exit 1
# else
  # echo "Recombination filtered alignment has already been bulit. Skipping this step..."
# fi

# if [ ! -f $outdir/filtered_core_aln.gubbins_filtered.snp_counts.txt ]; then
  # count_snps_in_core_alignment.py filtered_core_aln.gubbins_filtered.aln filtered_core_aln.gubbins_filtered.snp_counts.txt filtered_core_aln.gubbins_filtered.snp_counts.stats $threads || count_snps_in_core_alignment_fail=true
# else
  # echo "Core snp count matrix has already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/RAxML_bootstrap.filtered_core_aln.bootstrap ]; then
  # raxmlHPC-PTHREADS -T $threads -m GTRCAT -V -p 64855 -b 64855 -# 100 -s filtered_core_aln.filtered_polymorphic_sites.fasta -n filtered_core_aln.bootstrap || exit 1
# else
  # echo "Bootstrap RAxML trees have already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/RAxML_bipartitionsBranchLabels.filtered_core_aln.final ]; then
  # raxmlHPC-PTHREADS -T $threads -m GTRCAT -V -p 64855 -f b -t filtered_core_aln.final_tree.tre -z RAxML_bootstrap.filtered_core_aln.bootstrap -n filtered_core_aln.final || exit 1
# else
  # echo "Final RAxML tree has already been generated. Skipping this step..."
# fi

# if [ ! -f $outdir/RAxML_rootedTree.filtered_core_aln.finalrooted ]; then
  # raxmlHPC-PTHREADS -T $threads -m GTRCAT -V -f I -t RAxML_bipartitionsBranchLabels.filtered_core_aln.final -n filtered_core_aln.finalrooted || exit 1
# else
  # echo "Final RAxML rooted tree has already been generated. Skipping this step..."
# fi

# cd $outdir
# mkdir -p intermediate_files
# mv core.full.trimmed.aln \
# core.full.trimmed.aln.seq.joint.txt \
# core.branch_base_reconstruction.embl \
# core.filtered_polymorphic_sites.fasta \
# core.filtered_polymorphic_sites.phylip \
# core.final_tree.tre \
# core.per_branch_statistics.csv \
# core.gubbins_filtered.trimmed.aln \
# core.gubbins_filtered.trimmed.aln.reduced \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.0 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.1 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.2 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.3 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.4 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.5 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.0 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.1 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.2 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.3 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.4 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.5 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.6 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.6 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.7 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.8 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.9 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.10 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.11 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.12 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.7 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.8 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.9 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.10 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.11 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.12 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.13 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.13 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.14 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.15 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.16 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.17 \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.18 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.14 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.15 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.16 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.17 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.18 \
# RAxML_parsimonyTree.core.gubbins_filtered.trimmed.RUN.19 \
# RAxML_bestTree.core.gubbins_filtered.trimmed \
# RAxML_info.core.gubbins_filtered.trimmed \
# RAxML_log.core.gubbins_filtered.trimmed.RUN.19 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.0 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.1 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.2 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.3 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.4 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.5 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.6 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.7 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.8 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.9 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.10 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.11 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.12 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.13 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.14 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.15 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.16 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.17 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.18 \
# RAxML_result.core.gubbins_filtered.trimmed.RUN.19 \
# RAxML_bootstrap.core.gubbins_filtered.trimmed.bootstrap \
# RAxML_info.core.gubbins_filtered.trimmed.bootstrap \
# RAxML_bipartitions.core.gubbins_filtered.trimmed.final \
# RAxML_bipartitionsBranchLabels.core.gubbins_filtered.trimmed.final \
# RAxML_info.core.gubbins_filtered.trimmed.final \
# RAxML_info.core.gubbins_filtered.trimmed.finalrooted \
# intermediate_files

# if [ "$gubbins_drawer_fail" == true ]; then
  # echo "Gubbins drawer script failed. pdf of predicted recombinations for visualisation has not been generated"
# fi

# if [ "$count_snps_in_core_alignment_fail" == true ]; then
  # echo "Script to generate core snp count matrix has failed."
# fi

