#!/bin/bash
cd temp
javac -d classes -cp classes:opencsv-3.8.jar DrawSixFrame.java
for i in 1 2 3 4 5 6 9 10 11
do
	sixframe=$1.$i.sixframe
	lca=found_lcas_$1.$i.txt
	lineage=$1.$i.lineage
	java -cp .:classes:opencsv-3.8.jar drawsixframe.DrawSixFrame $i $lca $sixframe $2 > $1.$i.unknown.tex
	java -cp .:classes:opencsv-3.8.jar drawsixframe.DrawSixFrame $i $lca $sixframe $2 $lineage > $1.$i.known.tex  
done
