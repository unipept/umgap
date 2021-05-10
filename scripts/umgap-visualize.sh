#!/bin/sh
set -e
tmp="$(mktemp -d)"
trap "rm -rf '$tmp'" EXIT TERM

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

Usage: $0 -t [-r rank] <input>...
       $0 -w <input>
       $0 -u <input>...

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
while getopts c:r:wtu f; do
	case "$f" in
	c) configdir="$OPTARG" ;;
	r) rank="$OPTARG" ;;
	w) type="html" ;;
	t) type="csv" ;;
	u) type="url" ;;
	\?) crash "$USAGE" '' ;;
	esac
done
shift "$(( OPTIND - 1 ))"

[ "$#" -lt 1 ] && crash "$USAGE"
[ -z "$type" ] && crash "$USAGE"
[ "$type" = "html" -a "$#" -gt 1 ] && crash "$USAGE"

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
url)
	for file in "$@"; do
		filetype="$(file --mime-type "$file")" || crash "Could not determine filetype of '$file'."
		printf "%s: " "$file"
		if [ "$filetype" != "${filetype%gzip}" ]; then
			log "Inputfile '$file' is compressed"
			zcat "$file"
		else
			cat "$file"
		fi | umgap taxa2tree --url
	done
	;;
html)
	umgap taxa2tree < "$1" ;;
csv)
	inputfiles=""
	for file in "$@"; do
		filetype="$(file --mime-type "$file")" || crash "Could not determine filetype of '$file'."
		filename="$(printf '%s' "$file" | tr -c '[:alnum:].-' '_')"
		mkfifo "$tmp/$filename"
		if [ "$filetype" != "${filetype%gzip}" ]; then
			log "Inputfile '$file' is compressed"
			zcat "$file" > "$tmp/$filename" &
		else
			cat "$file" > "$tmp/$filename" &
		fi
		inputfiles="$inputfiles $tmp/$filename"
	done
	umgap taxa2freq -r "$rank" "$(getconfigdir)/$version/taxons.tsv" $inputfiles \
		| sed '1s_,[^,]*/_,_g'
	;;
esac


