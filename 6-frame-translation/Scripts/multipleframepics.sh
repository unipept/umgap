#!/bin/bash
mkdir temp_$1
cp ../DrawSixFrame/src/drawsixframe/DrawSixFrame.java temp_$1
for i in 1 2 3 4 5 6 9 10 11
do
	sixframe=$1.$i.sixframe
	lca=found_lcas_$1.$i.txt
	lineage=$1.$i.lineage
	cp ../SixFrameTransl/results/$2/$sixframe temp_$1
	cp ../SixFrameTransl/results/$2/$lca temp_$1
	cut -d "," -f 3 temp_$1/$lca | sed '1d' | while read line
	do 
		efetch -db taxonomy -id "$line" -format xml | xtract -pattern Taxon -element Lineage
	done > temp_$1/$lineage
done
