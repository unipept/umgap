#!/bin/sh
set -e

# print stuff to stderr and exits with fault
crash() {
	echo "$*" >&2
	exit 1
}

# ask a yes/no question.
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

if [ -z "$XDG_DATA_HOME" ]; then
	if [ -d "$HOME/Library/Application Support" ]; then
		data_default="$HOME/Library/Application Support/Unipept"
	elif [ -d "$HOME/.local/share" ]; then
		data_default="$HOME/.local/share/unipept"
	else
		data_default="$HOME/.unipept/data"
	fi
else
	data_default="$XDG_DATA_HOME/unipept"
fi


DATASERVER='https://unipept.ugent.be/system/umgap'
USAGE="
Setup the UMGAP by downloading the database and finding dependencies.

Usage: $0 [options]

Options:
  -c dir   The configuration directory. Defaults to '$config_default'.
  -f dir   Directory containing the FragGeneScan++ binary (FGSpp) and
           training data.
  -d dir   The data directory. The database files will be downloaded here.
           Defaults to '$data_default'.
  -y       Download all files of the latest version without asking.
"

while getopts c:f:d:y f; do
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
	configdir="$config_default"
	if [ -d "$configdir" ]; then
		echo "Using existing '$configdir' as config directory."
	else
		if ! yesno "Use '$configdir' as configuration directory?"; then
			crash "Please set a configuration directory with '-c'."
		fi
	fi
fi
if [ ! -d "$configdir" ]; then
	if mkdir -p "$configdir"; then
		echo "Created directory ${configdir}."
	else
		crash 'Could not create the configuration directory.'
	fi
fi

# If FGSpp is found and works, link it in the config directory.
if [ -z "$fgsppdir" ]; then
	if ! yesno 'FragGeneScan++ (-f) was not provided, only pipelines not using it will be available.'; then
		crash 'Pass the location of FragGeneScan++ with -f.'
	fi
else
	cd "$fgsppdir"
	if ! printf '>a\nA' | ./FGSpp -s stdin -o stdout -w 0 -t illumina_10 -c 240 > /dev/null; then
		crash 'Invoking FGSpp failed, follow installation instructions at https://github.com/unipept/FragGeneScanPlusPlus.'
	fi
	fgsppdir="$(pwd -P)" # make path absolute
	cd - > /dev/null
	if ln -f -s "$fgsppdir" "$configdir/FGSpp"; then
		echo 'Found, tested and remembered the FragGeneScan++ location.'
	else
		crash 'Could not link the FragGeneScan++ location.'
	fi
fi

# If no data directory is provided, use the default. Create it.
if [ -z "$datadir" ]; then
	datadir="$data_default"
	if [ -d "$datadir" ]; then
		echo "Using existing '$datadir' as data directory."
	else
		if ! yesno "Use '$datadir' as data directory?"; then
			crash "Please set a data directory with '-d'."
		fi
	fi
fi
if [ ! -d "$datadir" ]; then
	if mkdir -p "$datadir"; then
		echo "Created directory ${datadir}."
	else
		crash "Could not create the data directory."
	fi
fi
cd "$datadir"
datadir="$(pwd -P)" # make path absolute
cd - > /dev/null

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

if wget -V > /dev/null 2>&1; then
	get() { wget -q -O - "$1"; }
	size() { wget --spider --server-response "$1" 2>&1 | sed -n 's/ *[Cc]ontent-[Ll]ength: \([0-9]*\)\r*/\1/p' | human; }
	download() { wget -O "$2" "$1"; }
elif curl -V > /dev/null 2>&1; then
	get() { curl -s "$1"; }
	size() { curl -s -I "$1" | sed -n 's/ *[Cc]ontent-[Ll]ength: \([0-9]*\)\r*/\1/p' | human; }
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

echo
echo 'For any type of analysis, you need a taxonomony. For mapping tryptic or'
echo '9-mer peptides to it, you need the respective index file of the same'
echo 'version.'
echo

# Checking the files for this version.
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
	chmod 644 "$datadir/$version/taxons.tsv"
	ln -f -s "$datadir/$version/taxons.tsv" "$configdir/$version/taxons.tsv"
fi

if [ "$download_tryptic" = "true" ]; then
	download "$tryptic" "$datadir/$version/tryptic.fst"
	chmod 644 "$datadir/$version/tryptic.fst"
	ln -f -s "$datadir/$version/tryptic.fst" "$configdir/$version/tryptic.fst"
fi

if [ "$download_ninemer" = "true" ]; then
	download "$ninemer" "$datadir/$version/ninemer.fst"
	chmod 644 "$datadir/$version/ninemer.fst"
	ln -f -s "$datadir/$version/ninemer.fst" "$configdir/$version/ninemer.fst"
fi
