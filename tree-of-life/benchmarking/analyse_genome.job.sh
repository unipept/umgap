#!/bin/sh
#
#
#PBS -N normal
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=2:00:00
#PBS -m be
#

module load Ruby/2.2.1-intel-2014b
module load Python/3.3.2-ictce-4.1.13

# Export gempath
export PATH="$(ruby -rubygems -e 'puts Gem.user_dir')/bin:$PATH"

# make sure $PBS_ARRAYID is set, to avoid surprises
if [ -z "$PBS_ARRAYID" ]
then
  echo 'ERROR: $PBS_ARRAYID is not set, submit job(s) using "qsub -t <array expression>"'
  exit 1
fi

cd ..

# Do some parsing based on the arrayid
INPUTS_LIST_FILE=$HOME/metagenomics-scripts/tree-of-life/benchmarking/data/complete_assemblies.tsv
INPUT_LINE=$(sed -n "${PBS_ARRAYID}p" $INPUTS_LIST_FILE)
ASM_COL=$(echo $INPUT_LINE | awk -F\t '{print $9}')
ASM_ID=$(echo $NINTH_COL | sed 's/ .*//')

./benchmarking/analyse_genome.sh $ASM_ID -d $VSC_DATA -t $TMPDIR
