"""Calculates the LCA for a given protein"""

import sys
import subprocess
from itertools import zip_longest

from tree_of_life import get_tree, CLASSES

tree = get_tree()
unfound = set()
unmatched = set()
correct = 0
counter = 0


def get_lca(lineages):
    """Does the actual LCA calculation"""

    # Use -1 as fillvalue here, we'll filter it out later
    for i, taxons in enumerate(zip_longest(*lineages, fillvalue=-1)):
        # Remove the filling and the invalid taxons
        taxons_filtered = [t for t in taxons if t == -1 or t is None or t.valid_taxon]
        taxon_set = set(taxons_filtered) - set([-1])

        if CLASSES[i] == 'genus' or CLASSES[i] == 'species':
            taxon_set = taxon_set - set([None])

        if len(taxon_set) == 1:
            val = taxon_set.pop()
            if val:
                lca = val
        elif len(taxon_set) > 1:
            return lca

    return lca


def handle_lca(line):
    lineages = []
    lca = None

    print("Calculating LCA for {}".format(line))
    prot_result = subprocess.Popen("unipept pept2prot -s taxon_id {}".format(line), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for prot in prot_result.stdout.readlines()[1:]:
        taxon = tree.taxons[int(prot)]
        lineage = taxon.get_lineage(no_rank=False, invalid=False)
        lineages.append(lineage)

    for lineage in lineages:
        print([taxon.taxon_id if taxon else 0 for taxon in lineage])

    if not lineages:
        unfound.add(line)
        print("LCA: No LCA found for {}".format(line))
    else:
        lca = get_lca(lineages).taxon_id

        print("LCA: {}: {}".format(lca, tree.taxons[lca].name))

    return lca, lineages


def get_unipept_lca(line):
    unipept_lca = None

    result = subprocess.Popen("unipept pept2lca -s taxon_id -s taxon_name {}".format(line), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    unipept_result = result.stdout.readlines()
    if len(unipept_result) > 1:
        unipept_lca = unipept_result[1].decode('utf-8').split(',')
        print("Unipept LCA: {}: {}".format(unipept_lca[0], unipept_lca[1]))
    else:
        print("Unipept LCA: None LCA found for {}".format(line))

    return unipept_lca


for line in sys.stdin:
    lca, lineages = handle_lca(line)
    unipept_lca = get_unipept_lca(line)

    if not lineages and not unipept_lca:
        correct = correct + 1
    elif lineages and unipept_lca and lca == int(unipept_lca[0]):
        correct = correct + 1
    else:
        unmatched.add(line)

    counter = counter + 1
    print()
    print("=================== {}% ==================".format(counter/3983*100))
    print("Unmatched: {}, {}".format(len(unmatched), ', '.join(unmatched)))
    print("Correct: {}, total: {}, accuracy: {}%".format(correct, counter, correct/counter*100))
    print("=================== {}% ==================".format(counter/3983*100))
    print()

print("Unfound: {}, {}".format(len(unfound), ', '.join(unfound)))
