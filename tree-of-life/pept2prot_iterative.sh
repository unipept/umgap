#!/bin/bash

diff results/all.txt results/processed.txt  --suppress-common-lines  | cut -d' ' -f2 | tail -n+2 > results/unprocessed.txt

i=1
while read file
do
  echo "Processing file $i $file"
  cat "results/fastafiles/$file" | unipept pept2prot --host scruffy.ugent.be --extra > "results/pept2prot_fastafiles/$file"
  i=$((i+1))
done < results/unprocessed.txt
