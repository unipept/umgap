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

# Create dirs
dir="/tmp/$ass_id/"

rm -rf $dir
mkdir $dir


#  get the complete sequence and process it with:
#     - prot2pept
#     - peptfilter
python ./entrez/ass2pept.py $ass_id > "$dir/peptides.fst"

# get the proteins uniprot ids which occur in the genome
python ./entrez/ass2seqacc.py $ass_id | entrez/seqacc2protid.sh > "$dir/uniprot_protein_ids.txt"

# analyse the complete sequence with:
#     - pept2lca2lca
#     - pept2prot2filter2lca

# check wether resulting taxons come from the correct lineage and reports how many do

# spit out some statistics about the found lcas


