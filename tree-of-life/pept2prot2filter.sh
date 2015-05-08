#!/bin/bash

usage(){
  echo "Usage: cat pept2prot.fst | $0 uniprot_identifiers_file" >&2
  exit 1
}

(($# != 1)) && usage

grep -v -f $1
