#!/bin/bash


FILES=results/pept2prot_fastafiles/*

i=1
for file in $FILES
do
  filename=$(basename $file)

  echo "Processing file $i $filename"
  grep -v "^>|\([^,]*\),.*,[^,]*\1[^,]*$" $file > "/work/tnnaesse/$filename"
  i=$((i+1))
done
