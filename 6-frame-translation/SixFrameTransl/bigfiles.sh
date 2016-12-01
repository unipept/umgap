#!/bin/bash
sixframe=$1.sixframe
cat $sixframe | prot2pept | peptfilter > temp
lines=$(wc -l temp | sed 's/^[ ]*//' | cut -d " " -f 1)
for i in $(seq 500 500 8062018) 
do
	head -n $i temp | tail -n 500 | unipept pept2lca >> $1.lcas.txt
done
