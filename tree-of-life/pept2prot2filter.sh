#!/bin/bash

usage(){
  echo "Usage: $0 uniprot_identifiers_file pept2prot_fasta_file"
  exit 1
}

(($# != 2)) && usage

grep -v -f $1 $2
