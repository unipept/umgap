#!/bin/sh
set -e
self="$$"
tmp="$(mktemp -d)"

# Preparing exit
children=""

trapf() {
	rm -rf "$tmp"
	[ -n "$children" ] && kill $children
}

trap trapf EXIT KILL

# What to do if the XDG standard isn't there...
if [ -z "$XDG_CONFIG_HOME" ]; then
	if [ -d "$HOME/Library/Preferences/" ]; then
		config_default="$HOME/Library/Preferences/Unipept"
	elif [ -d "$HOME/.config" ]; then
		# why weren't the XDG variables set then?
		config_default="$HOME/.config/unipept"
	else
		config_default="$HOME/.unipept"
	fi
else
	config_default="$XDG_CONFIG_HOME/unipept"
fi

USAGE="
Some default configurations for the UMGAP.

Usage: $0 -1 <fasta> [-t <type>] [-z] -o <output>
       $0 -1 <fastq1> -2 <fastq2> [-t <type>] [-z] -o <output>

Where:
  <fasta>   A (optionally GZIP-compressed) FASTA file of reads.
  <fastq1>  A (optionally GZIP-compressed) FASTQ file of reads (first end).
  <fastq2>  A (optionally GZIP-compressed) FASTQ file of reads (second end).
  <type>    The type of analysis, any of:
              - max-sensitivity
              - high-sensitivity
              - tryptic-sensitivity
              - tryptic-precision
              - high-precision (the default)
              - max-precision
  <output>  The output file. Use '-' to write to stdout.

Options:
  -z        GZIP-compress before writing the following output file.
  -c dir    The configuration directory. Defaults to '$config_default'.

Note: you can call this script with multiple samples, by just repeating above
arguments. It will then share files loaded in memory in between samples to save
time:

  $0 -1 <fastaA> -o <outputA> -1 <fastq1B> -2 <fastq2B> -t <typeB> -o <outputB> ...
"

# =========================================================================== #
#  Some functions.
# =========================================================================== #

# Logging
log() {
	[ -z "$VERBOSE" ] && return
	printf "log: %s\n" "$*" >&2
}

debug() {
	[ -z "$VERBOSE" -a -z "$DEBUG" ] && return
	printf "debug: %s\n" "$*" >&2
}

# print stuff to stderr and exits with fault
crash() {
	debug "encountered error"
	echo "$*" >&2
	finish 1
}

# check if a filename contains tabs or newlines
checkname() {
	debug "checking if filename '$1' is valid"
	[ "${1}" = "${1%	*}" ] || crash "Tabs in filenames are unsupported."
}

# function to fetch the configuration directory
configdir=""
getconfigdir() {
	if [ -n "$configdir" ]; then
		echo "$configdir"
	elif [ -d "$config_default" ]; then
		echo "$config_default"
	elif [ -d /etc/umgap ]; then
		echo /etc/umgap
	else
		crash "No configuration directory found. Please run umgap-setup or use the '-c' argument."
	fi
}

# =========================================================================== #
#  Environmental checks.
# =========================================================================== #

debug "checking if umgap is installed"
if ! umgap -V > /dev/null; then
	crash 'Cannot find the umgap executable. Please ensure it is installed and located in your $PATH.'
fi

# =========================================================================== #
#  Argument parsing.
# =========================================================================== #

debug "parsing the arguments"

count=0
fgspp_used="false"
tryptic_used="false"
ninemer_used="false"
samples=""

# defaults
compress="false"
type="high-precision"
tmpfiles=""
while getopts c:1:2:t:zo: f; do
	case "$f" in
	c) configdir="$OPTARG" ;;
	z) compress="true" ;;
	1) checkname "$OPTARG"; infile1="$OPTARG" ;;
	2) checkname "$OPTARG"; infile2="$OPTARG" ;;
	t) checkname "$OPTARG"; type="$OPTARG" ;;
	o) checkname "$OPTARG"; outfile="$OPTARG" ;;
	\?) crash "$USAGE" ;;
	esac

	if [ "$f" = "o" ]; then # o is the last option of each series
		count="$(( count + 1 ))"
		log "Checking sample $count:"

		case "$type" in
		"max-sensitivity") ninemer_used="true" ;;
		"high-sensitivity") ninemer_used="true" ;;
		"tryptic-sensitivity") fgspp_used="true" ; tryptic_used="true" ;;
		"tryptic-precision") fgspp_used="true" ; tryptic_used="true" ;;
		"high-precision") fgspp_used="true" ; ninemer_used="true" ;;
		"max-precision") fgspp_used="true" ; ninemer_used="true" ;;
		*) crash "Unrecognized type: '$type'."
		esac
		log "- Using '$type'."

		if [ -n "$infile1" ]; then
			log "- found input file 1: '$infile1'"
			filetype="$(file --mime-type "$infile1")" || \
				crash "Could not determine filetype of '$infile1'."
			if [ "$filetype" != "${filetype%gzip}" ]; then
				log "- input file 1 is compressed"
				mkfifo "$tmp/gunzip1-$count"
				{ zcat "$infile1" > "$tmp/gunzip1-$count" || kill "$self"; } > /dev/null &
				infile1="$tmp/gunzip1-$count"
			fi
		fi

		if [ -n "$infile2" ]; then
			log "- found input file 2: '$infile2'"
			filetype="$(file --mime-type "$infile2")" || \
				crash "Could not determine filetype of '$infile2'."
			if [ "$filetype" != "${filetype%gzip}" ]; then
				log "- input file 2 is compressed"
				mkfifo  "$tmp/gunzip2-$count"
				{ zcat "$infile2" > "$tmp/gunzip2-$count" || kill "$self"; } > /dev/null &
				infile2="$tmp/gunzip2-$count"
			fi
		fi

		if [ -n "$infile2" ]; then
			if [ -n "$infile1" ]; then
				log "- found two input files, assuming paired-end FASTQ"
				mkfifo "$tmp/infile-$count"
				{ umgap fastq2fasta "$infile1" "$infile2" > "$tmp/infile-$count" || kill "$self"; } > /dev/null &
				infile="$tmp/infile-$count"
			else
				crash "Encountered a second input file without a first."
			fi
		elif [ -n "$infile1" ]; then
			infile="$infile1"
		else
			crash "Encountered an output file without input files."
		fi

		if [ "$outfile" = "-" ]; then
			mkfifo "$tmp/outfile-$count"
			cat "$tmp/outfile-$count" &
			outfile="$tmp/outfile-$count"
		fi

		if [ "$compress" = "true" ]; then
			mkfifo "$tmp/gzip-$count"
			{ gzip - < "$tmp/gzip-$count" > "$outfile" || kill "$self"; } > /dev/null &
			outfile="$tmp/gzip-$count"
		fi

		samples="$(printf '%s%s\t%s\t%s\t' "$samples" "$type" "$infile" "$outfile")"

		log ""
		
		# reset for next
		compress="false"
		infile1=
		infile2=
		type="high-precision"
		outfile=
	fi
done

if [ "$count" -eq 0 ]; then
	echo "$USAGE"
	finish
fi

# =========================================================================== #
#  Environmental checks.
# =========================================================================== #

if [ "$fgspp_used" = "true" ]; then
	debug "checking if FGSpp is installed"
	if ! [ -h "$(getconfigdir)/FGSpp" ]; then
		crash "FragGeneScan++ not found. Please run umgap-setup."
	fi
fi

versions="$(find -H "$(getconfigdir)" -mindepth 1 -maxdepth 1 \
                 -printf '%P\n' | sort -n)"
for candidate in $versions; do
	[ ! -h "$(getconfigdir)/$candidate/taxons.tsv" ] && continue
	[ "$tryptic_used" = "true" -a ! -h "$(getconfigdir)/$candidate/tryptic.fst" ] && continue
	[ "$ninemer_used" = "true" -a ! -h "$(getconfigdir)/$candidate/ninemer.fst" ] && continue
	version="$candidate"
done
[ -n "$version" ] || crash "No data version found valid for all samples. Please run umgap-setup."
debug "using version '$version'"

# =========================================================================== #
#  The actual pipeline
# =========================================================================== #

fgspp() {
	"$(getconfigdir)/FGSpp/FGSpp" -s stdin -o stdout -w 0 \
	    -r "$(getconfigdir)/FGSpp/train" -t "illumina_10" -p 4 -c 2
}

taxons="$(getconfigdir)/$version/taxons.tsv"
ninemers="$(getconfigdir)/$version/ninemer.fst"
tryptics="$(getconfigdir)/$version/tryptic.fst"

if [ "$ninemer_used" = "true" ]; then
	log "Loading index in memory."
	socket="$tmp/socket"
	{ umgap prot2kmer2lca -m -o -s "$socket" "$ninemers" || kill "$self"; } > /dev/null &
	ninemer_pid="$!"
	children="$(pgrep -P "$!") $ninemer_pid"
	while [ ! -S "$socket" ] && sleep 1; do true; done
fi

# for each sample
count=0
while [ -n "${samples}" ]; do
	type="${samples%%	*}" ; samples="${samples#*	}"
	infile="${samples%%	*}" ; samples="${samples#*	}"
	outfile="${samples%%	*}" ; samples="${samples#*	}"

	count="$(( count + 1 ))"
	log "Running sample ${count}."

	case "$type" in
	"max-sensitivity")
		umgap translate -a                             | # translate 6 frames
		socat - UNIX-CONNECT:"$socket"                 | # map to taxa
		umgap seedextend -g1 -s2                       | # seedextend filter
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l1 -m rmq -a mrtl "$taxons"   ;; # aggregate
	"high-sensitivity")
		umgap translate -a                             | # translate 6 frames
		socat - UNIX-CONNECT:"$socket"                 | # map to taxa
		umgap seedextend -g1 -s3                       | # seedextend filter
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l1 -a hybrid -f 0.25 "$taxons";; # aggregate
	"tryptic-sensitivity")
		fgspp                                          | # gene prediction
		umgap prot2tryp2lca    -l9 -L45 "$tryptics"    | # map to taxa
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l1 -m rmq -a mrtl "$taxons"   ;; # aggregate
	"tryptic-precision")
		fgspp                                          | # gene prediction
		umgap prot2tryp2lca    -l9 -L45 "$tryptics"    | # map to taxa
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l5 -m rmq -a mrtl "$taxons"   ;; # aggregate
	"high-precision")
		fgspp                                          | # gene prediction
		socat - UNIX-CONNECT:"$socket"                 | # map to taxa
		umgap seedextend -g1 -s3                       | # seedextend filter
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l2 -a lca\* "$taxons"         ;; # aggregate
	"max-precision")
		fgspp                                          | # gene prediction
		socat - UNIX-CONNECT:"$socket"                 | # map to taxa
		umgap seedextend -g1 -s4                       | # seedextend filter
		umgap uniq -d /                                | # join paired ends
		umgap taxa2agg -l5 -a lca\* "$taxons"         ;; # aggregate
	esac < "$infile" > "$outfile"
done
