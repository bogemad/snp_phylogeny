#!/bin/bash

base_path=$PWD
reads_dir=`readlink -f raw_data/reads`

reference=`find raw_data/reference_sequence -type f`
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
  hpc=$2
fi
id_threshold=$3

echo ""
echo "Reads directory = $reads_dir"
echo "Reference sequence file = $reference"
echo "Output directory = $outdir"
echo "Process threads = $threads"
echo "Minimum reference coverage threshold = $id_threshold%"
echo ""

ref_files=`find raw_data/reference_sequence -type f | wc -l`
if [[ $ref_files > 1 ]]; then
echo "Make sure there is only one reference file in raw_data/reference_sequence directory."
exit 1
fi

if [[ "$reference" == *.gbk || "$reference" == *.gb || "$reference" == *.gbff ]]; then
  refsuff=".gbk"
elif [[ "$reference" == *.fasta || "$reference" == *.fa || "$reference" == *.fna ]]; then
  refsuff=".fasta"
else
  echo "Can't id reference sequence format. Ensure you have a Genbank (.gb, .gbk, .gbff) or Fasta (.fasta,.fa,.fna) formatted sequence."
  exit 1
fi

if [ ! -f $outdir/reference$refsuff ]; then
  remove_degenerate_bases.py $reference $outdir/reference $threads || exit 1
else
  echo "Reference sequence has already been processed. Skipping this step..."
fi

reference=$outdir/reference$refsuff
cd $reads_dir

# for file in `find $reads_dir -type f`; do
  # name=`get_samplenames.py $file`
  # if [ $name == 'fail' ]; then echo "Failed to determine samplename for file path: $name. Check if this is a fastq or fastq.gz file."; exit 1; fi
  # if [ ! -f $outdir/$name/$name.aligned.fa ]; then
    # snippy --prefix $name --cpus $threads --outdir $outdir/$name --ref $reference --peil $file || (echo "Alignment failed! Check if $filename is correct fastq file." && exit 1)
  # else
    # echo "$name reads have already been aligned. Skipping this step..."
  # fi
# done
exit_on_end=false
#if [ "$hpc" == "true" ]
#then
#  for file in `find -L $reads_dir -type f`
#  do
#    name=`get_samplenames.py $file`
#    if [ ! -f $outdir/$name/$name.aligned.fa ]; then
#      exit_on_end=true
#      if [ "${#name}" -gt 15 ]
#      then
#        short_name=${name:0:15}
#      else
#        short_name=$name
#      fi
#      sed "s~xxxshortjobnamexxx~$short_name~g" $base_path/scripts/qsub_multisnippy.sh | sed "s~xxxbasepathxxx~$base_path~g" | sed "s~xxxfilexxx~$file~g" | sed "s~xxxoutdirxxx~$outdir~g" | sed "s~xxxreferencexxx~$reference~g" > $temp/qsub_$name.sh
#      qsub $temp/qsub_$name.sh
#      rm -rf $temp/qsub_$name.sh
#    fi
#  done
#else
  parallel -j $threads --progress multisnippy.sh {} $outdir $reference ::: `find -L $reads_dir -type f`
#fi

if [ "$exit_on_end" == "true" ]
then
echo "Snippy jobs have been submitted to the HPC smallq for processing. Please restart the pipeline when these jobs have been completed."
exit
fi

cd $outdir
rm $reference


if [ ! -f $outdir/core.full.aln ]; then
  snippy-core * || exit 1
  mkdir -p ../excluded_sequences/poor_ref_alignment
  ditch_distant_core_sequences.py $id_threshold
  snippy-core * || exit 1
else
  echo "Core alignment has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/core.full.trimmed.aln ]; then
  replace_Ns_with_gaps.py core.full.aln $threads core.full.trimmed.aln || exit 1
else
  echo "Core alignment has already been trimmed. Skipping this step..."
fi

if [ ! -f $outdir/total.snp_counts.txt ]; then
  count_snps_in_core_alignment.py core.full.trimmed.aln total.snp_counts.txt $threads
else
  echo "Total snp count matrix has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/core.final_tree.tre ]; then
  echo "Running recombination filter..."
  run_gubbins.py -i 10 -f 15 -p core -c $threads core.full.trimmed.aln || exit 1
else
  echo "Recombination filtering by gubbins has already been done. Skipping this step..."
fi

if [ ! -f $outdir/core.recombination_predictions.pdf ]; then
  gubbins_drawer.py -o core.recombination_predictions.pdf -t core.final_tree.tre core.recombination_predictions.embl || gubbins_drawer_fail=true
else
  echo "Recombination prediction visualisation has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/core.gubbins_filtered.aln ]; then
  filter_recombinations.py core.recombination_predictions.gff core.full.trimmed.aln core.gubbins_filtered.aln || exit 1
else
  echo "Recombination filtered alignment has already been bulit. Skipping this step..."
fi

if [ ! -f $outdir/core.gubbins_filtered.snp_counts.txt ]; then
  count_snps_in_core_alignment.py core.gubbins_filtered.aln core.gubbins_filtered.snp_counts.txt $threads || count_snps_in_core_alignment_fail=true
else
  echo "Core snp count matrix has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/core.gubbins_filtered.trimmed.aln ]; then
  remove_gaps_Ns.py core.gubbins_filtered.aln $threads core.gubbins_filtered.trimmed.aln
else
  echo "Gap filtered alignment has already been bulit. Skipping this step..."
fi


if [ ! -f $outdir/RAxML_bestTree.core.gubbins_filtered.trimmed ]; then

  raxmlHPC-PTHREADS -T $threads -m GTRGAMMA -p 64855 -# 20 -s core.gubbins_filtered.trimmed.aln -n core.gubbins_filtered.trimmed || exit 1
else
  echo "Initial RAxML tree has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/RAxML_bootstrap.core.gubbins_filtered.trimmed.bootstrap ]; then
  raxmlHPC-PTHREADS -T $threads -m GTRGAMMA -p 64855 -b 64855 -# 100 -s core.gubbins_filtered.trimmed.aln -n core.gubbins_filtered.trimmed.bootstrap || exit 1
else
  echo "Bootstrap RAxML trees have already been generated. Skipping this step..."
fi

if [ ! -f $outdir/RAxML_bipartitionsBranchLabels.core.gubbins_filtered.trimmed.final ]; then
  raxmlHPC-PTHREADS -T $threads -m GTRGAMMA -p 64855 -f b -t RAxML_bestTree.core.gubbins_filtered.trimmed -z RAxML_bootstrap.core.gubbins_filtered.trimmed.bootstrap -n core.gubbins_filtered.trimmed.final || exit 1
else
  echo "Final RAxML tree has already been generated. Skipping this step..."
fi

if [ ! -f $outdir/RAxML_rootedTree.core.gubbins_filtered.trimmed.finalrooted ]; then
  raxmlHPC-PTHREADS -T $threads -m GTRGAMMA -f I -t RAxML_bipartitionsBranchLabels.core.gubbins_filtered.trimmed.final -n core.gubbins_filtered.trimmed.finalrooted || exit 1
else
  echo "Final RAxML rooted tree has already been generated. Skipping this step..."
fi

if [ "gubbins_drawer_fail" == true ]; then
  echo "Gubbins drawer script failed. pdf of predicted recombinations for visualisation has not been generated"
fi

if [ "count_snps_in_core_alignment_fail" == true ]; then
  echo "Script to generate core snp count matrix has failed."
fi

