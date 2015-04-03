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
tmp_dir="/tmp"

# Create dirs
dir="$tmp_dir/$ass_id"

#rm -rf $dir
mkdir -p $dir

# get the taxon ID of the assembly
tax_id=$(python ./entrez/ass2taxid.py $ass_id)

#  get the complete sequence and process it with:
#     - prot2pept
#     - peptfilter
if [ ! -f "$dir/peptides.fst" ]; then
  python ./entrez/ass2pept.py $ass_id > "$dir/peptides.fst"
fi

# get the proteins uniprot ids which occur in the genome
if [ ! -f "$dir/uniprot_protein_ids.txt" ]; then
  python ./entrez/ass2seqacc.py $ass_id | entrez/seqacc2protid.sh > "$dir/uniprot_protein_ids.txt"
fi

# analyse the complete sequence with and
# check wether resulting taxons come from the correct lineage
#     - pept2lca2lca
unipept pept2lca -i "$dir/peptides.fst" | python3 pept2lca2lca.py -c $tax_id

#     - pept2prot2filter2lca
unipept pept2prot -i "$dir/peptides.fst"| ./pept2prot2filter.sh "$dir/uniprot_protein_ids.txt" | python3 pept2prot2filter2lca.py -c $tax_id


# spit out some statistics about the found lcas
