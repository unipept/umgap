#!/bin/sh
set -e
tmp="$(mktemp -d)"
trap "rm -rf '$tmp'" EXIT KILL

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
Visualizing data with the UMGAP.

Usage: $0 -i <input> [-r rank] -t
       $0 -i <input> -w
       $0 -i <input> -u

Where:
  <input>   A (optionally GZIP-compressed) FASTA file of taxa.

Options:
  -t        Output a CSV frequency table on species rank.
  -w        Output an HTML webpage of an interactive visualization.
  -u        Print a shareable URL to a online interactive visualisation.
  -c dir    The configuration directory. Defaults to '$config_default'.
  -r rank   Set the rank for the CSV frequency table (default: species).
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
	exit 1
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
#  Argument parsing.
# =========================================================================== #

debug "parsing the arguments"

rank="species"
while getopts i:c:r:wtu f; do
	case "$f" in
	c) configdir="$OPTARG" ;;
	i) inputfile="$OPTARG" ;;
	r) rank="$OPTARG" ;;
	w) type="html" ;;
	t) type="csv" ;;
	u) type="url" ;;
	\?) crash "$USAGE" '' ;;
	esac
done

[ -z "$inputfile" ] && crash "$USAGE"
[ -z "$type" ] && crash "$USAGE"

filetype="$(file --mime-type "$inputfile")" || \
	crash "Could not determine filetype of '$inputfile'."
if [ "$filetype" != "${filetype%gzip}" ]; then
	log "Inputfile is compressed"
	mkfifo "$tmp/gunzip"
	zcat "$inputfile" > "$tmp/gunzip" &
	inputfile="$tmp/gunzip"
fi

# =========================================================================== #
#  Environmental checks.
# =========================================================================== #

debug "checking if umgap is installed"
if ! umgap -V > /dev/null; then
	crash 'Cannot find the umgap executable. Please ensure it is installed and located in your $PATH.'
fi

debug "checking if we have a taxons file for the frequency table"
if [ "$type" = "csv" ]; then
	versions="$(find -H "$(getconfigdir)" -mindepth 1 -maxdepth 1 \
	                 -printf '%P\n' | sort -n)"
	for candidate in $versions; do
		[ ! -h "$(getconfigdir)/$candidate/taxons.tsv" ] && continue
		version="$candidate"
	done
	[ -n "$version" ] || crash "No taxon table found for frequency counting. Please run umgap-setup."
	debug "using version '$version'"
fi

# =========================================================================== #
#  The actual visualization code
# =========================================================================== #

case "$type" in
url) umgap taxa2tree < "$inputfile" --url ;;
html) umgap taxa2tree < "$inputfile" ;;
csv) grep -v '^>' "$inputfile" | umgap taxa2freq -r "$rank" "$(getconfigdir)/$version/taxons.tsv" ;;
esac
