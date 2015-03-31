#!/bin/bash

usage(){
  echo "Usage: cat pept2prot | $0 uniprot_identifiers_file"
  exit 1
}

(($# != 2)) && usage

grep -v -f $1
