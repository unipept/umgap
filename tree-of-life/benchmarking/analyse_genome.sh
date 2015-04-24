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
  echo "Usage: $0 [refseq assembly id]"
  exit 1
}


(($# != 1)) && usage

ass_id=$1
dir=$(mktemp -d -t $ass_id)

#rm -rf $dir
mkdir -p $dir

# get the taxon ID of the assembly
tax_id=$(python3 ./entrez/ass2taxid.py $ass_id)

#  get the complete sequence and process it with:
#     - prot2pept
#     - peptfilter
python3 ./entrez/ass2pept.py $ass_id > "$dir/peptides.fst"

# get the proteins uniprot ids which occur in the genome
python3 ./entrez/ass2seqacc.py $ass_id | entrez/seqacc2protid.sh > "$dir/uniprot_protein_ids.txt"

# analyse the complete sequence with and
# check wether resulting taxons come from the correct lineage
#     - pept2lca2lca
unipept pept2lca -i "$dir/peptides.fst" | python3 pept2lca2lca.py -c $tax_id > "$dir/pept2lca2lca.fst"

#     - pept2prot2filter2lca
#unipept pept2prot -i "$dir/peptides.fst"| ./pept2prot2filter.sh "$dir/uniprot_protein_ids.txt" | python3 pept2prot2filter2lca.py -c $tax_id > "$dir/pept2prot2filter2lca.fst"


# spit out some statistics about the found lcas
