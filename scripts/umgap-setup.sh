#!/bin/sh
# - install rust? or ask to
# - install umgap
# - ask configdir (or use XDG/unipept or ~/.unipept)
# - install FGSpp (or ask to) and find it (symlink in configdir)
# - download taxons (ask location? to XDG_DATA_HOME/unipept) (symlink in configdir) (VERSION!)
# - download indices (ask first: if, where) (symlink in configdir) (VERSION!)

# print stuff to stderr and exits with fault
crash() {
	echo "$*" >&2
	exit 1
}

DATASERVER='https://unipept.ugent.be/system/umgap'
USAGE="
Setup the UMGAP by downloading the database and finding dependencies.

Usage: $0 [options]

Options:
  -c dir   The configuration directory. Defaults to '$XDG_CONFIG_HOME/unipept'.
  -f dir   Directory containing the FragGeneScan++ binary (FGSpp) and
           training data.
  -d dir   The data directory. The database files will be downloaded here.
           Defaults to '$XDG_DATA_HOME/unipept'.
  -y       Download all files of the latest version without asking.
"

while getopts c:f:t:n:y f; do
	case "$f" in
	c) configdir="$OPTARG" ;;
	f) fgsppdir="$OPTARG" ;;
	d) datadir="$OPTARG" ;;
	y) download_all="true" ;;
	\?) crash "$USAGE" ;;
	esac
done

# If no config directory is provided, use the default. Create it.
if [ -z "$configdir" ]; then
	echo "Config directory not set, using $XDG_CONFIG_HOME/unipept."
	configdir="$XDG_CONFIG_HOME/unipept"
fi
if mkdir -p "$configdir"; then
	echo "Created directory ${configdir}."
else
	crash 'Could not create the configuration directory.'
fi

# If FGSpp is found and works, link it in the config directory.
if [ -z "$fgsppdir" ]; then
	echo 'FragGeneScan++ was not provided, only pipelines not using it will be available.'
else
	cd "$fgsppdir"
	if ! printf '>a\nA' | ./FGSpp -s stdin -o stdout -w 0 -t illumina_10 -m 1024 > /dev/null; then
		crash 'Invoking FGSpp failed, follow installation instructions at https://github.com/unipept/FragGeneScanPlusPlus.'
	fi
	fgsppdir="$PWD" # make path absolute
	cd - > /dev/null
	if ln -f -s "$fgsppdir" "$configdir/FGSpp"; then
		echo 'Found, tested and remembered the FragGeneScan++ location.'
	else
		crash 'Could not link the FragGeneScan++ location.'
	fi
fi

# If no data directory is provided, use the default. Create it.
if [ -z "$datadir" ]; then
	echo "Data directory not set, using $XDG_DATA_HOME/unipept."
	datadir="$XDG_DATA_HOME/unipept"
fi
if mkdir -p "$datadir"; then
	echo "Created directory ${datadir}."
	cd "$datadir"
	datadir="$PWD" # make path absolute
	cd - > /dev/null
else
	crash "Could not create the data directory."
fi

# How can we download stuff
human() {
	if numfmt --version > /dev/null; then
		numfmt --to=iec-i --suffix=B
		return
	fi

	number="$(cat -)"
	for unit in "B" "KiB" "MiB" "GiB"; do
		if [ "$number" -lt 9999 ]; then
			printf '%d%s' "$number" "$unit"
			return
		fi
		number="$(( number / 1024 ))"
	done
	printf '%d%s' "$number" "$unit"
}

if wget -V > /dev/null; then
	get() { wget -q -O - "$1"; }
	size() { wget --spider --server-response "$1" 2>&1 | sed -n 's/ *[Cc]ontent-[Ll]ength: \([0-9]*\)/\1/p' | human; }
	download() { wget -O "$2" "$1"; }
elif curl -V > /dev/null; then
	get() { curl -s "$1"; }
	size() { curl -I "$1" | sed -n 's/ *[Cc]ontent-[Ll]ength: \([0-9]*\)/\1/p' | human; }
	download() { curl -o "$2" "$1"; }
else
	crash 'Neither curl nor wget is available to download data.'
fi

# Retrieving the latest version from unipept
echo 'Checking the latest version on the server.'
version="$(get "$DATASERVER/latest")"
if [ "$?" != 0 ]; then
	crash 'Could not retrieve version from server.'
fi
echo "Latest version is ${version}."

# Checking the files for this version.
yesno() {
	if [ "$download_all" = "true" ]; then
		true
	else
		question="$1 [y]/n "
		read -p "$question" answer
		while [ "$answer" != "y" -a "$answer" != "" -a "$answer" != "n" ]; do
			echo "$answer is not a valid option."
			read -p "$question" answer
		done
		[ "$answer" != "n" ]
	fi
}

taxons="$DATASERVER/$version/taxons.tsv"
if [ -h "$configdir/$version/taxons.tsv" ]; then
	echo "Taxonomy is available at the latest version."
elif yesno "Would you like to download the taxonomy from ${version} ($(size "$taxons"))?"; then
	download_taxons="true"
fi

tryptic="$DATASERVER/$version/tryptic.fst"
if [ -h "$configdir/$version/tryptic.fst" ]; then
	echo "Tryptic index is available at the latest version."
elif yesno "Would you like to download the tryptic index from ${version} ($(size "$tryptic"))?"; then
	download_tryptic="true"
fi

ninemer="$DATASERVER/$version/ninemer.fst"
if [ -h "$configdir/$version/ninemer.fst" ]; then
	echo "9-mer index is available at the latest version."
elif yesno "Would you like to download the 9-mer index from ${version} ($(size "$ninemer"))?"; then
	download_ninemer="true"
fi

if [ "$download_taxons" = "true" -o "$download_ninemer" = "true" -o "$download_ninemer" = "true" ]; then
	mkdir -p "$datadir/$version"
	mkdir -p "$configdir/$version"
fi

if [ "$download_taxons" = "true" ]; then
	download "$taxons" "$datadir/$version/taxons.tsv"
	ln -s "$datadir/$version/taxons.tsv" "$configdir/$version/taxons.tsv"
fi

if [ "$download_tryptic" = "true" ]; then
	download "$tryptic" "$datadir/$version/tryptic.fst"
	ln -s "$datadir/$version/tryptic.fst" "$configdir/$version/tryptic.fst"
fi

if [ "$download_ninemer" = "true" ]; then
	download "$ninemer" "$datadir/$version/ninemer.fst"
	ln -s "$datadir/$version/ninemer.fst" "$configdir/$version/ninemer.fst"
fi
