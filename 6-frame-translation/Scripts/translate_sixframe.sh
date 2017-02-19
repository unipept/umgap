#!/bin/bash
# A shell script to use the java script to generate a sixframe translation of a given fasta file with a nucleotide sequence
# This script assumes the java script is located in a directory called SixFrameTransl with the same parent directory as the folder where this script is located

usage(){ echo "Syntax: $(basename $0) [-t tt1-tt2-...] [-k k:kmer-index] [-p] FASTAfile outputPrefix outputLocation" 1>&2; exit 1; }

translationTables="1 2 3 4 5 6 9 10 11"
while getopts ":t:k:p" opt; do
	case $opt in
    		t) translationTables=$(echo "$OPTARG" | tr '-' ' ') 	;;
		k) kmer=true
		   k=$(echo "$OPTARG"|cut -d':' -f1) 
		   kmer_index=$(echo "$OPTARG"| cut -d':' -f2)	;;
		p) tripept=true ;;
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
outputLocation=$3

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


# outputLocation=$(cd $outputLocation | pwd)


cd "$(dirname $0)/../SixFrameTransl"

for i in $translationTables
do
	mvn exec:java -Dexec.mainClass="com.mycompany.sixframetransl.translate" -Dexec.args="$FASTAfile $i False" | egrep -v '(INFO|WARNING)'  > "$outputLocation/$outputPrefix.$i.sixframe"
	if [ "$tripept" ]
	then
		cat "$outputLocation/$outputPrefix.$i.sixframe" | prot2pept | peptfilter |unipept pept2lca > "$outputLocation/found_lcas_$outputPrefix.$i.txt"
	fi
	if [ "$kmer" ]
	then 
		if [[ -z "$(which prot2kmer)" ]]; then
			echo 'prot2kmer is not in your PATH' >&2
        		exit 2
		fi

		if [[ -z "$(which pept2lca)" ]]; then
			echo 'pept2lca is not in your PATH' >&2
			exit 2
		fi
		cat "$outputLocation/$outputPrefix.$i.sixframe" | prot2kmer -k $k | fasta-uniq -kw | pept2lca -i "$kmer_index"  > "$outputLocation/$outputPrefix.$i.${k}mer.lca"
	fi
done

