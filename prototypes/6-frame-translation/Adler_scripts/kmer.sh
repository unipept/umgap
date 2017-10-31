#!/bin/bash

export PATH="/data/felix/metagenomics/repo/unipept/target/release/:$PATH"

name="$1"
sixframe="$1.sixframe"

if [[ -z "$(which prot2kmer)" ]]; then
	echo 'prot2kmer is not in your PATH' >&2
        exit 2
fi

if [[ -z "$(which pept2lca)" ]]; then
	echo 'pept2lca is not in your PATH' >&2
	exit 2
fi

kmer_index="/tmp/tmpfs/9mers.index"
taxon="$(dirname $0)/taxons.9mer.tsv"

RUST_BACKTRACE=1
#cat "$sixframe" | prot2kmer -k 9 | fasta-uniq -kw | pept2lca -i -t $taxon "$kmer_index"  > "$name.9mer.lca"
#cat "$sixframe" | prot2kmer -k 9 | fasta-uniq -kw | pept2lca "$kmer_index"  > "$name.9mer.lca"
cat "$sixframe" | prot2kmer -k 9 | fasta-uniq -kw | pept2lca "$kmer_index" > /dev/null
