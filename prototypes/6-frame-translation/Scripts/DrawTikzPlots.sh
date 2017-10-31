#!/bin/bash
# Script to run the Java-code that makes the plots

usage(){ echo "Syntax: $(basename $0) [-k k(Integer)] [-t tt1-tt2-...] [-l organismName] outputLocation  name  inputLocation  title  [proteinPositions]" 1>&2; exit 1;}

translationTables="1 2 3 4 5 6 9 10 11"
while getopts ":t:k:l:" opt; do
        case $opt in
                t) translationTables=$(echo "$OPTARG" | tr '-' ' ')     ;;
		k) k=$OPTARG ;;
		l) organism=$OPTARG ;;
                :) usage ;;
                \?) usage
  esac
done
shift $((OPTIND - 1))

if [ $# -lt 4 -o $# -gt 5 ]
then
	usage
fi

name="$2"
title="$4"
if [ $# -eq 5 ]
then
	proteinPos="$5"
fi

if [ ! -d "$1" -o ! -r "$1" ]
then
	echo "directory $1 not found" 1>&2
	exit 2
else
	outputLocation=$(cd $1 ;pwd)
	if [ ! -d "$outputLocation/$name" ]
	then
		mkdir "$outputLocation/$name"
	fi
fi

if [ ! -d "$3" -o ! -r "$3" ]
then 
	echo "directory $3 not found" 1>&2
        exit 2
else
	inputLocation=$(cd $3 ;pwd)
fi

echo "\documentclass[10pt,a4paper,landscape]{report}
\usepackage{amsfonts}
\usepackage{graphicx}
\usepackage[dvipsnames]{xcolor}
\usepackage{tikz}
\usepackage[left=2cm,right=2cm,top=3cm,bottom=3cm]{geometry}
\usetikzlibrary{decorations.pathreplacing}
\begin{document}" 


for i in $translationTables
do
	sixframe=$name.$i.sixframe
	if [ "$k" ]
	then
		lca=$name.$i.${k}mer.lca
	else
		#lca=found_lcas_$name.$i.txt
		lca=$name.$i.pept.lca
	fi
	if [ "$organism" ]
	then
		if [ "$k" ]
		then
			lineage=$inputLocation/$name.$i.kmer.lineage
		else
			lineage=$inputLocation/$name.$i.pept.lineage
		fi
		if [ ! -f "$lineage" ]
		then
			if [ "$k" ]
			then
				egrep -v '^>' $inputLocation/$lca | cut -d ',' -f2 | while read line
	                        do
                                	if [ ! $line -eq 1 ]
					then
						unipept taxonomy -a -s "*_name" $line | tail -1 | cut -d',' -f2- | tr -s ',' |sed 's/,$//'
					fi 
        	                done > "$lineage"
			else
				#cut -d "," -f 3 $inputLocation/$lca | sed '1d' | while read line
				egrep -v '^>' $inputLocation/$lca | cut -d ',' -f2 | while read line
				do 
					if [ ! $line -eq 1 ]
                                        then
                                                unipept taxonomy -a -s "*_name" $line | tail -1 | cut -d',' -f2- | tr -s ',' |sed 's/,$//'
                                        fi
				done > "$lineage"
			fi
		fi
		organism_lineage="$(esearch -db taxonomy -query "$organism" | efetch -format xml | xtract -pattern Taxon -element Lineage); $organism"
	fi	
done

cd ~/Documents/Unief/Thesis/Thesis/ScoreReads
rm -r classes
mkdir classes
javac -d classes -cp classes:../opencsv-3.8.jar:taxonomy_score.txt $(find src -name *.java)

for i in $translationTables
do
	sixframe=$inputLocation/$name.$i.sixframe
	if [ "$k" ]
        then
                lca=$inputLocation/$name.$i.${k}mer.lca
		tex=$name.$i.${k}mer
        else
                #lca=$inputLocation/found_lcas_$name.$i.txt
        	lca=$inputLocation/$name.$i.pept.lca
		tex=$name.$i.pept
	fi
	if [ "$organism" ]
	then
		if [ "$k" ]
                then
                        lineage=$inputLocation/$name.$i.kmer.lineage
                else
                        lineage=$inputLocation/$name.$i.pept.lineage
                fi
		java -cp .:classes:../opencsv-3.8.jar:taxonomy_score.txt scorereads.ScoreReads ${k:-0} "$title, translation table $i" $sixframe $lca ${proteinPos:-:} $lineage "${organism_lineage}" > $outputLocation/$name/$tex.known.tex
		echo "\input{$tex.known.tex}" 
		echo
	else
		java -cp .:classes:../opencsv-3.8.jar:taxonomy_score.txt scorereads.ScoreReads ${k:-0} "$title, translation table $i" $sixframe $lca ${proteinPos:-:} > $outputLocation/$name/$tex.unknown.tex 
		echo "\input{$tex.unknown.tex}" 
		echo	
fi
done
echo "\end{document}"
exit 0
