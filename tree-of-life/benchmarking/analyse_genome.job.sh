#!/bin/sh
#
#
#PBS -l walltime=2:00:00
#PBS -m n
#

module load Ruby/2.2.1-intel-2015a
module load Python/3.4.3-intel-2015a
module load Biopython/1.65-intel-2015a-Python-3.4.3

ASM_ID=${asm_id}

cd $HOME/unipept-metagenomics-scripts/tree-of-life
./benchmarking/analyse_genome.sh $ASM_ID -d $VSC_DATA -t $VSC_DATA -r "$VSC_SCRATCH/.rmqdatadir"
