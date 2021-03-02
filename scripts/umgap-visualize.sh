#!/bin/sh
set -e

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
fi

USAGE="
Visualizing data with the UMGAP.

Usage: $0 -i <input> -o <output>
       $0 -i <input> -u

Where:
  <input>   A (optionally GZIP-compressed FASTA file of taxa.
  <output>  The HTML output file. Use '-' to write to stdout.

Options:
  -u        Instead of writing to a file, print a URL to an
            online visualisation
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
	printf "log: %s\n" "$*" >&2
}

# cleans up temporary files and exits
tmpfiles=""
finish() {
	debug "removing temporary files"
	while [ -n "${tmpfiles%%*}" ]; do
		debug "removing ${tmpfiles%%*}"
		rm -f "${tmpfiles%%*}"
		tmpfiles="${tmpfiles#*}"
	done
	debug "quitting"
	exit "${1:-0}"
}

getfifo() {
	name="$(gettempname)"
	mkfifo "$name"
	printf '%s' "$name"
}

gettempname() {
	name="$(mktemp)"
	rm "$name"
	printf '%s' "$name"
}

# print stuff to stderr and exits with fault
crash() {
	debug "encountered error"
	echo "$*" >&2
	finish 1
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

while getopts i:o:u f; do
	case "$f" in
	i) inputfile="$OPTARG" ;;
	o) outputfile="$OPTARG" ;;
	u) url="yes" ;;
	\?) crash "$USAGE" ''
	esac
done

[ -z "$inputfile" ] && crash "$USAGE"
[ -n "$outputfile" -a -n "$url" ] && crash "$USAGE"
[ -z "$outputfile" -a -z "$url" ] && crash "$USAGE"

filetype="$(file --mime-type "$inputfile")" || \
	crash "Could not determine filetype of '$inputfile'."
if [ "$filetype" != "${filetype%gzip}" ]; then
	log "Inputfile is compressed"
	fifo="$(getfifo)"
	zcat "$inputfile" > "$fifo" &
	inputfile="$fifo"
fi

if [ "$outputfile" = "-" ]; then
	fifo="$(getfifo)"
	cat "$fifo" &
	outfile="$fifo"
fi

# =========================================================================== #
#  The actual visualization code
# =========================================================================== #

if [ -n "$url" ]; then
    umgap taxa2tree < "$inputfile" --url
else
    umgap taxa2tree < "$inputfile" > "$outputfile"
fi
