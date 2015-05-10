#!/bin/bash
# usage: cat peptides.fst | commonpeptfilter.sh pept2lca.fst

grep -v -F <(grep '\(no rank\|superkingdom\)' $1 | cut -f2 -d, | sort | uniq)
