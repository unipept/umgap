#!/bin/sh
#PBS -N build-index
#PBS -q long
#PBS -m abe
#PBS -l nodes=1:ppn=8
#PBS -l vmem=200gb

# Loading the required modules
module load Cargo

# Check the current directory
pushd "$PBS_O_WORKDIR"

# Running the build-index
./target/release/build-index <(zcat "$VSC_DATA_VO_USER/data/tables") "$VSC_DATA_VO_USER/data/fst/index.bin"

# Reset the directory
popd


