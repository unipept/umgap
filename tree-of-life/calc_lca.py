#!/usr/bin/env python

import sys
import subprocess

from tree_of_life import get_tree

def get_lca(lineages):
    index = 1
    last = lineages[0][0]

    while index < len(lineages[0]):
        first = lineages[0][index]
        for lineage in lineages:
            if first != lineage[index]:
                return last

        last = first
        index = index + 1

    return last


tree = get_tree()
unfound = set()

for line in sys.stdin:

    p = subprocess.Popen("unipept pept2prot -s taxon_id {}".format(line), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    lineages = []

    for line in p.stdout.readlines()[1:]:
        taxon = tree.taxons[int(line)]
        lineages.append(taxon.get_lineage())

    for lineage in lineages:
        print(lineage)

    if not lineages:
        unfound.add(line)
    else:
        if len(lineages) == 1:
            lca = lineages[0][-1]
        else:
            lca = get_lca(lineages)

        print(lca)
        print(tree.taxons[lca].name)


print("Unfound: {}".format(unfound))
