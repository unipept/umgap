#!/bin/bash
#PBS -N build-index
#PBS -q long
#PBS -m abe
#PBS -l nodes=1:ppn=3
#PBS -l vmem=200gb

# Running the build-index
"$PBS_O_WORKDIR/build-index" <(zcat "$VSC_DATA_VO_USER/data/tables/sequences.tsv.gz") "$VSC_DATA_VO_USER/data/fst/index.bin"

