#!/bin/bash

# Analyses one genome:

#  gets the assemblies references of the genome
#  gets the complete sequence
#  gets the proteins uniprot ids which occur in the genome
#  processes the complete sequence with:
#      - prot2pept
#      - peptfilter
#  analyses the complete sequence with:
#      - pept2lca2lca
#      - pept2prot2filter2lca
#  checks wether resulting taxons come from the correct lineage and reports how many do
#  spits out some statistics about the found lcas

set -eu

usage() {
  echo "Usage: $0 [refseq assembly id] [-d datadir] [-t tempdir] [-r rmqdatadir]" >&2
  exit 1
}


(($# < 1)) && usage

# Save directory of the analysis script to know where to find the others
dir="$(dirname "$0")"

asm_id=$1 && shift

# Create a tmpdir and a datadir
tmpdir=$(mktemp -d -t "$asm_id.XXXXXXXXXX")
datadir=$tmpdir
rmqdatadir=""

while getopts "d:t:r:" opt
do
  case $opt in
    d)
      datadir="$OPTARG/$asm_id"
      ;;
    t)
      tmpdir="$OPTARG/$asm_id"
      ;;
    r)
      rmqdatadir="-r $OPTARG"
      ;;
    ?)
      usage
      ;;
  esac
done

mkdir -p $datadir
mkdir -p $tmpdir

echo "Writing data to $datadir"
echo "Writing tempdata to $tmpdir"

# get the taxon ID of the assembly
tax_id=$(python3 $dir/../entrez/asm2taxid.py $asm_id)

#  get the complete sequence and process it with:
#     - prot2pept
#     - peptfilter

echo "Getting peptides"
if [ ! -f "$tmpdir/peptides.fst" ]
then
  python3 $dir/../entrez/asm2pept.py $asm_id > "$tmpdir/peptides.fst"
fi

# analyse the complete sequence with and
# check wether resulting taxons come from the correct lineage
#     - pept2lca2lca
echo "Executing pept2lca2lca"
unipept pept2lca -i "$tmpdir/peptides.fst" \
  | tee "$datadir/pept2lca.fst" \
  | python3 $dir/../pept2lca2lca.py -c $tax_id $rmqdatadir > "$datadir/pept2lca2lca.fst"


echo "Getting uniprot ids"
# get the proteins uniprot ids which occur in the genome
if [[ -s "$tmpdir/uniprot_protein_ids.txt" ]]
then
  python3 $dir/../entrez/asm2seqacc.py $asm_id | python3 $dir/../entrez/seqacc2protid.py > "$tmpdir/uniprot_protein_ids.txt"

  # Check if the file isn't empty
  if [[ -s "$tmpdir/uniprot_protein_ids.txt" ]]
  then
    echo "ERROR: It seems that the uniprot_protein_ids.txt file is empty. Pept2prot2filter has no use anymore now." >&2
    exit 1
  fi
fi

#     - pept2prot2filter2lca
echo "Executing pept2prot2filter2lca"
unipept pept2prot -i "$tmpdir/peptides.fst" \
  | tee "$datadir/pept2prot.fst" \
  | $dir/../pept2prot2filter.sh "$tmpdir/uniprot_protein_ids.txt" \
  | tee "$datadir/pept2prot2filter.fst" \
  | python3 $dir/../pept2prot2filter2lca.py -c $tax_id $rmqdatadir > "$datadir/pept2prot2filter2lca.fst"


# spit out some statistics about the found lcas


#rm -rf $tmpdir
