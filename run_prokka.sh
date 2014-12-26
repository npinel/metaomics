#!/bin/bash -l

module load BLAST
module load HMMER
module load Infernal
module load parallel/20131022-goolf-1.4.10

date

/home/users/snarayanasamy/prokka-1.7/bin/prokka-1.7.1\
  --force\
  --outdir $2
  --cpus 6\
  --metagenome\
  --norrna\
  $1

date

