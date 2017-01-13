#!/bin/bash
# A shell script to use the java script to generate a sixframe translation of a given fasta file with a nucleotide sequence
# This script assumes the java script is located in a directory called SixFrameTransl with the same parent directory as the folder where this script is located

usage(){ echo "Syntax: $(basename $0) [-t tt1-tt2-...] FASTAfile outputPrefix outputLocation" 1>&2; exit 1; }

translationTables="1 2 3 4 5 6 9 10 11"
while getopts ":t:" opt; do
	case $opt in
    		t) translationTables=$(echo "$OPTARG" | tr '-' ' ') 	;;
		:) usage ;;
    		\?) usage
  esac
done
shift $((OPTIND - 1))

if [ $# -ne 3 ]
then 
	usage
fi

FASTAfile="$1"
outputPrefix="$2"
outputLocation="$3"

if [ ! -r "$FASTAfile" -o ! -f "$FASTAfile" ]
then
	echo "$FASTAfile can\'t be found or read " 1>&2
	exit 2
fi

if [ ! -d "$outputLocation" ]
then
	echo "location $outputLocation not found" 1>&2
	exit 2
fi


if [[ "$outputLocation" = "." ]]
then 
	outputLocation="$(pwd)"
fi


cd "$(dirname $0)/../SixFrameTransl"

for i in "$translationTables"
do
	mvn exec:java -Dexec.mainClass="com.mycompany.sixframetransl.translate" -Dexec.args="$FASTAfile $i False" | egrep -v '(INFO|WARNING)'  > $outputPrefix.$i.sixframe
	cat $outputPrefix.$i.sixframe | prot2pept | peptfilter |unipept pept2lca > found_lcas_$outputPrefix.$i.txt
	mv $outputPrefix.$i.sixframe $outputLocation
	mv found_lcas_$outputPrefix.$i.txt $outputLocation
done

