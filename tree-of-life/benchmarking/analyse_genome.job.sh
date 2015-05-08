#!/bin/sh
#
#
#PBS -N debug
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=0:59:58
#PBS -m be
#

module load Ruby/2.2.1-intel-2015a
module load Python/3.4.3-intel-2015a
module load Biopython/1.65-intel-2015a-Python-3.4.3

# make sure $PBS_ARRAYID is set, to avoid surprises
if [ -z "$PBS_ARRAYID" ]
then
  echo 'ERROR: $PBS_ARRAYID is not set, submit job(s) using "qsub -t <array expression>"'
  exit 1
fi


# Do some parsing based on the arrayid
INPUTS_LIST_FILE=$HOME/unipept-metagenomics-scripts/tree-of-life/benchmarking/data/complete_assemblies.tsv
ASM_ID=$(sed -n "${PBS_ARRAYID}p" $INPUTS_LIST_FILE | awk -F'\t' '{print $9}' | sed 's/ .*//')

cd $HOME/unipept-metagenomics-scripts/tree-of-life
./benchmarking/analyse_genome.sh $ASM_ID -d $VSC_DATA -t $TMPDIR -r "$VSC_SCRATCH/.rmqdatadir"
