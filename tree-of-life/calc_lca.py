#!/usr/bin/env python

import sys
import subprocess
from itertools import zip_longest

from tree_of_life import get_tree

def get_lca(lineages):

    lca = 1 # We always start with the root
    for x in zip_longest(*lineages, fillvalue=0):
        a = set(x) - set([0])
        if len(a) == 1:
            lca = a.pop()
        elif len(a) > 1:
            return lca

    return lca


tree = get_tree()
unfound = set()
correct = 0
counter = 0

for line in sys.stdin:
    lineages = []

    print("Calculating LCA for {}".format(line))
    prot_result = subprocess.Popen("unipept pept2prot -s taxon_id {}".format(line), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for prot in prot_result.stdout.readlines()[1:]:
        taxon = tree.taxons[int(prot)]
        lineage = taxon.get_lineage(no_rank=False, invalid=False)
        lineages.append(lineage)

    for lineage in lineages:
        print(lineage)

    if not lineages:
        unfound.add(line)
        print("LCA: No LCA found for {}".format(line))
    else:
        lca = get_lca(lineages)

        print("LCA: {}: {}".format(lca, tree.taxons[lca].name))

    unipept_result = subprocess.Popen("unipept pept2lca -s taxon_id -s taxon_name {}".format(line), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = unipept_result.stdout.readlines()
    if len(result) > 1:
        unipept_lca = result[1].decode('utf-8').split(',')
        print("Unipept LCA: {}: {}".format(unipept_lca[0], unipept_lca[1]))
    else:
        print("Unipept LCA: None LCA found for {}".format(line))

    if not lineages and not result:
        correct = correct + 1
    elif lineages and result and lca == int(unipept_lca[0]):
        correct = correct + 1
    counter = counter + 1
    print()
    print("=================== {}% ==================".format(counter/3983*100))
    print()

print("Unfound: {}, {}".format(len(unfound), ', '.join(unfound)))
print("Correct: {}, total: {}, accuracy: {}%".format(correct, counter, correct/counter*100))
