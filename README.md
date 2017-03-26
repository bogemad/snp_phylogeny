# snp_phylogeny
## A base for running snp based phylogenetic analysis

##### Download SRA reads from an exported enterobase table (reads will be deposited into the raw_data/reads directory):

`./download_enterobase_SRA_reads.sh <path_to_enterobase_table> <number_of_concurrent_downloads_to_request>`

Add input reads and reference sequence to the raw_data/reads and raw_data/reference_sequence directories. Only one reference sequence can be used currently.

##### To run snp_phylogeny pipeline:

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
  threads                     number of cpus requested from hpc and number of processes used by pipeline (as 
                              above) (usually equiv. to number of cpu cores).
  memory                      RAM requested from the hpc - use the format Xgb where X is the about of RAM 
                              to be requested.
  coverage_cutoff_threshold   % coverage of reference sequence (0-100%) used to reject a sample. Samples 
                              lower than this threshold will be excluded from phylogenetic pipeline steps.
  ```
