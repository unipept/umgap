""" Given a file of pept2lca entries, gets the LCA for every fasta_header.

Expects a list grouped first by the FASTA header and then by peptide for speed reasons
Not necessarily ordered, just grouped is fine
"""

import sys
from itertools import groupby
from operator import itemgetter

import argparse

from lca_calculators.tree_based_lca import Tree_LCA_Calculator

# Input argument parsing
parser = argparse.ArgumentParser(description='Calculate the LCAs for a given fastafile after pept2lca')
parser.add_argument('-c', '--check-against', dest='reference_taxon_id', type=int, help='check against a taxon_id')
args = parser.parse_args()


calculator = Tree_LCA_Calculator()

print("fasta_header,taxon_id,taxon_name,taxon_rank", end="")
if args.reference_taxon_id:
    print(",on_lineage")
else:
    print()


def reduce_per_fasta(lines):
    lca = calculator.calc_lca([int(line[2]) for line in lines], allow_no_rank=False)
    taxon = calculator.tree.taxons[lca]

    return taxon


def parse(input_):
    splits = (line.rstrip().split(',') for line in input_)
    fasta_groups = ((k, g) for k, g in groupby(splits, itemgetter(0)))
    reduced_fastas = ((k, reduce_per_fasta(group)) for k, group in fasta_groups)

    for fasta_header, taxon in reduced_fastas:

        print("{},{},{},{}".format(fasta_header, taxon.taxon_id, taxon.name, taxon.rank), end="")

        # We want to check against a reference taxon
        if args.reference_taxon_id:
            on_lineage = calculator.calc_lca([taxon.taxon_id, args.reference_taxon_id], allow_no_rank=True) == args.reference_taxon_id
            print(",{}".format(int(on_lineage)))
        else:
            print()


# Skip the header
next(sys.stdin)
parse(sys.stdin)

calculator.cleanup()
