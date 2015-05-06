#!/bin/sh
#
#
#PBS -N normal
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=2:00:00
#PBS -m be
#

module load Ruby
module load Python

# make sure $PBS_ARRAYID is set, to avoid surprises
if [ -z "$PBS_ARRAYID" ]
then
  echo 'ERROR: $PBS_ARRAYID is not set, submit job(s) using "qsub -t <array expression>"'
  exit 1
fi

cd ..

INPUTS_LIST_FILE=$HOME/metagenomics-scripts/tree-of-life/benchmarking/data/complete_assemblies.tsv
INPUT_LINE=`sed -n "${PBS_ARRAYID}p" $INPUTS_LIST_FILE`

./benchmarking/analyse_genome.sh $1 -d $VSC_DATA -t $TMPDIR
