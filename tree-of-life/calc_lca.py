#!/usr/bin/env python

import sys

from tree_of_life import get_tree

def get_lca(lineages):
    index = 1
    last = lineages[0][0]

    while True:
        first = lineages[0][index]
        for lineage in lineages:
            if first != lineage[index]:
                return last

        last = first
        index = index + 1


tree = get_tree()
lineages = {}

next(sys.stdin) # skip the header
for line in sys.stdin:
    taxon =  tree.taxons[int(line)]
    lineages[taxon] = taxon.get_lineage()

for lineage in lineages:
    print(lineages.get(lineage))

lca = get_lca(list(lineages.values()))
print(lca)
print(tree.taxons[lca].name)
