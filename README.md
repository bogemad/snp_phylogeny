# snp_phylogeny
## A base for running snp based phylogenetic analysis

### Installation

Download and install with:

```
git clone https://github.com/bogemad/snp_phylogeny.git <dirname>
cd <dirname>
make
```

where:

```
dirname = anything!!! Whatever you wish to call your analysis!
```

The pipeline will install all dependancies required for the anaysis.

### Running the pipeline

During installation raw_data/reference_sequence, raw_data/reads, analysis_results and logs directories will be generated. When run (using run.sh or run-hpc.sh, see below), snp_phylogeny will search for a reference sequence and reads in the raw_data/reference_sequence and raw_data/reads directories. Each fasta (assembled genome) or fastq (sequencing reads) file will be mapped to the reference sequence using snippy. Following this, alignments are concatenated using snippy-core and total snp counts calculated. Recombination events are consequently filtered with gubbins, core snp counts calculated and maximum likelihood trees built (and bootstapped) with RAxML.

To run the pipeline, add your reads and reference sequence to the raw_data/reads and raw_data/reference_sequence directories. Only one reference sequence can be used currently. 

Results will be deposited into the analysis_results directory.

##### Download SRA reads from an exported enterobase table (reads will be deposited into the raw_data/reads directory):

`./download_enterobase_SRA_reads.sh <path_to_enterobase_table> <number_of_concurrent_downloads_to_request>`

##### To run snp_phylogeny pipeline:

With reads and reference sequence in the raw_data directory enter:

`./run.sh <threads> <coverage_cutoff_threshold>`

where:

```
  threads                     number of processes used by the pipeline (usually equiv. to number of cpu cores).
  coverage_cutoff_threshold   % coverage of reference sequence (0-100%) used to reject a sample. Samples lower 
                              than this threshold will be excluded from phylogenetic pipeline steps.
  ```
  
##### To run on a PBS-type job submission system (like those use on some hpc's):

`./run-hpc.sh <threads> <memory> <coverage_cutoff_threshold>`

where:

```
  threads                     number of cpus to be requested from hpc and number of processes used by pipeline (as 
                              above).
  memory                      RAM requested from the hpc - use the format Xgb where X is the about of RAM 
                              to be requested.
  coverage_cutoff_threshold   % coverage of reference sequence (0-100%) used to reject a sample. Samples 
                              lower than this threshold will be excluded from phylogenetic pipeline steps.
  ```

###### Cleanup scripts

Run `./cleanup.sh` or `./cleanup-hpc.sh` scripts to remove all core alignment, snp count and RAxML files and return exluded sequences to the reads and analysis_results directories (Especially useful while you are tweaking the coverage threshold to suit your data).

### Results

Once the pipeline is completed results files can be found in the analysis_results directory. A directory for each reads file is generated and includes files generated by snippy (see https://github.com/tseemann/snippy for further details). Data from these alignments are used to generate the core genome alignments. Files for core genome alignments are found at the base of the analysis results directory and include:

| File | Description | Generated by | See for details |
| --- | --- | --- | --- |
| core.aln | A core SNP alignment in fasta format | snippy | https://github.com/tseemann/snippy |
| core.full.aln | A whole genome SNP alignment (includes invariant sites) | snippy | https://github.com/tseemann/snippy |
| core.tab | Tab-separated columnar list of core SNP sites with alleles and annotations | snippy | https://github.com/tseemann/snippy |
| core.txt | Tab-separated columnar list of alignment/core-size statistics | snippy | https://github.com/tseemann/snippy |
| core.vcf | a vcf formatted file of of core SNP sites | snippy | https://github.com/tseemann/snippy |
| total.snp_counts.txt | tab-separated table of total SNPs between pairwise combinations of all samples | snp_phylogeny | Right here! |
| core.recombination_predictions.embl | Recombination predictions in EMBL tab file format | gubbins | https://github.com/sanger-pathogens/gubbins |
| core.recombination_predictions.gff | Recombination predictions in GFF format | gubbins | https://github.com/sanger-pathogens/gubbins |
| core.summary_of_snp_distribution.vcf | VCF file summarising the distribution of SNPs | gubbins | https://github.com/sanger-pathogens/gubbins |
| core.recombination_predictions.pdf | Visual graphic of recombination predictions | gubbins | https://github.com/sanger-pathogens/gubbins |
| core.gubbins_filtered.aln | core genome alignment with recombination filtered by gubbins | snp_phylogeny | Right here! |
| core.gubbins_filtered.snp_counts.txt | tab-separated table of total SNPs between pairwise combinations of all samples | snp_phylogeny | Right here! |
| RAxML_rootedTree.core.gubbins_filtered.trimmed.finalrooted | final ultrametric rooted phylogenetic tree | RAxML | http://sco.h-its.org/exelixis/web/software/raxml/index.html |

Intermediate files generated during the pipeline are stored in the intermediate_files directory withing analysis_results. Explanations for these are coming soon...

