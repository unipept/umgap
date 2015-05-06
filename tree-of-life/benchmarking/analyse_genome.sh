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


usage() {
  echo "Usage: $0 [refseq assembly id] [-d datadir] [-t tempdir]"
  exit 1
}


(($# < 1)) && usage

asm_id=$1

tmpdir=$(mktemp -d -t "$asm_id")
datadir=$tmpdir

if [ "$2" == "-d" ]
then
  datadir="$3/$asm_id"
fi

if [ "$4" == "-t" ]
then
  tmpdir="$5/$asm_id"
fi

mkdir -p $datadir
mkdir -p $tmpdir

echo "Writing data to $datadir"
echo "Writing tempdata to $tmpdir"

exit

# get the taxon ID of the assembly
tax_id=$(python3 ./entrez/asm2taxid.py $asm_id)

#  get the complete sequence and process it with:
#     - prot2pept
#     - peptfilter
python3 ./entrez/asm2pept.py $asm_id > "$tmpdir/peptides.fst"

# get the proteins uniprot ids which occur in the genome
python3 ./entrez/asm2seqacc.py $asm_id | entrez/seqacc2protid.sh > "$datadir/uniprot_protein_ids.txt"

# analyse the complete sequence with and
# check wether resulting taxons come from the correct lineage
#     - pept2lca2lca
unipept pept2lca -i "$tmpdir/peptides.fst" | tee "$datadir/pept2lca.fst" | python3 pept2lca2lca.py -c $tax_id > "$datadir/pept2lca2lca.fst"

#     - pept2prot2filter2lca
#unipept pept2prot -i "$tmpdir/peptides.fst"| ./pept2prot2filter.sh "$tmpdir/uniprot_protein_ids.txt" | python3 pept2prot2filter2lca.py -c $tax_id > "$datadir/pept2prot2filter2lca.fst"


# spit out some statistics about the found lcas


#rm -rf $tmpdir
