#!/bin/bash
mkdir temp
cp ../DrawSixFrame/src/drawsixframe/DrawSixFrame.java temp
for i in 1 2 3 4 5 6 9 10 11
do
	sixframe=$1.$i.sixframe
	lca=found_lcas_$1.$i.txt
	lineage=$1.$i.lineage
	cp ../DrawSixFrame/$sixframe temp
	cp ../DrawSixFrame/$lca temp
	cut -d "," -f 3 temp/$lca | sed '1d' | while read line
	do 
		efetch -db taxonomy -id "$line" -format xml | xtract -pattern Taxon -element Lineage
	done > temp/$lineage
done
