#!/bin/sh
#
#
#PBS -N normal
#PBS -o output.file
#PBS -e error.file
#PBS -l walltime=2:00:00
#PBS -m be
#

./benchmarking/analyse_genome.sh $1
