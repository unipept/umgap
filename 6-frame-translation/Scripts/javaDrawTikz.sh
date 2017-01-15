#!/bin/bash
cp ../DrawSixFrame/src/drawsixframe/DrawSixFrame.java temp_$1
cp opencsv-3.8.jar temp_$1
cd temp_$1
echo "\documentclass[10pt,a4paper,landscape]{report}
\usepackage[latin1]{inputenc}
\usepackage[dutch]{babel}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{tikz}
\usepackage[left=2cm,right=2cm,top=3cm,bottom=3cm]{geometry}
\author{Aranka Steyaert}
\begin{document}"
mkdir classes
javac -d classes -cp classes:opencsv-3.8.jar DrawSixFrame.java
for i in 1 2 3 4 5 6 9 10 11
do
	sixframe=$1.$i.sixframe
	lca=found_lcas_$1.$i.txt
	lineage=$1.$i.lineage
	esearch -db taxonomy -query "$3" | efetch -format xml | xtract -pattern Taxon -element Lineage > organism.lineage
	java -cp .:classes:opencsv-3.8.jar drawsixframe.DrawSixFrame "$3" organism.lineage $i $lca $sixframe $2 > $1.$i.unknown.tex
	echo "\include{$1.$i.unknown}"
	java -cp .:classes:opencsv-3.8.jar drawsixframe.DrawSixFrame "$3" organism.lineage $i $lca $sixframe $2 $lineage > $1.$i.known.tex  
	echo "\include{$1.$i.known}"
done
echo "\end{document}"
