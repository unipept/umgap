"""Calculates the LCA for a given protein"""

import sys
import subprocess
from itertools import zip_longest

from tree_of_life import get_tree, CLASSES

TREE = get_tree()
MATCHES = []
UNFOUND = set()
UNMATCHED = set()
CORRECT = 0
COUNTER = 0


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
    """Gets the lineages from the tree and performs the LCA calculation"""
    lineages = []
    lca = None

    # Get the prots from Unipept
    prot_result = subprocess.Popen(
        "unipept pept2prot -s taxon_id {}".format(line),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    # Get the lineage for each prot
    for prot in prot_result.stdout.readlines()[1:]:
        taxon = TREE.taxons[int(prot)]
        lineage = taxon.get_lineage(no_rank=False, invalid=False)
        lineages.append(lineage)

    for lineage in lineages:
        print([taxon.taxon_id if taxon else 0 for taxon in lineage])

    # If we have a result, get the LCA, otherwise, add it to the unfound ones
    if lineages:
        lca = get_lca(lineages).taxon_id
        print("LCA: {}: {}".format(lca, TREE.taxons[lca].name))
    else:
        UNFOUND.add(line)
        print("LCA: No LCA found for {}".format(line))

    return lca, lineages


def get_unipept_lca(line):
    """Gets the Unipept LCA"""

    unipept_lca = None

    # Get the LCA from Unipept
    result = subprocess.Popen(
        "unipept pept2lca -s taxon_id -s taxon_name {}".format(line),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    unipept_result = result.stdout.readlines()

    # Hande the result
    if len(unipept_result) > 1:
        unipept_lca = unipept_result[1].decode('utf-8').split(',')
        print("Unipept LCA: {}: {}".format(unipept_lca[0], unipept_lca[1]))
    else:
        print("Unipept LCA: None LCA found for {}".format(line))

    return unipept_lca


def print_lcas(filename):
    lineages = [TREE.taxons[prot].get_lineage(no_rank=False, invalid=False) for prot in MATCHES if prot]

    # Remove all the children from the taxons
    for lineage in lineages:
        for taxon in lineage:
            if taxon:
                taxon.children = set()
                taxon.count = 0
                taxon.self_count = 0

    # Add the children from the found lineages back
    for taxons in zip_longest(*lineages, fillvalue=-1):
        taxons_set = set(taxons) - set([None]) - set([-1])

        # Add the correct children
        for taxon in taxons_set:
            valid_parent = taxon.get_valid_parent()
            if taxon.parent_id != taxon.taxon_id:
                valid_parent.children.add(taxon)

            # Also set the self_count
            taxon.self_count = MATCHES.count(taxon.taxon_id)

    # Percolate the counts up the tree
    lineages[0][0].add_counts()

    # Print the result
    with open(filename, "wb") as f:
        lineages[0][0].to_json(f)


for line in sys.stdin:
    print("Calculating LCA for {}".format(line))

    lca, lineages = handle_lca(line)
    unipept_lca = get_unipept_lca(line)

    if (not lineages and not unipept_lca) or \
            (lineages and unipept_lca and lca == int(unipept_lca[0])):
        CORRECT += 1
        MATCHES.append(lca)
    else:
        UNMATCHED.add(line)

    COUNTER += 1
    print()
    print("=================== {}% ==================".format(COUNTER/3983*100))
    print("Unmatched: {}, {}".format(len(UNMATCHED), ', '.join(UNMATCHED)))
    print("Correct: {}, total: {}, accuracy: {}%".format(CORRECT, COUNTER, CORRECT/COUNTER*100))
    print("=================== {}% ==================".format(COUNTER/3983*100))
    print()

print("Unfound: {}, {}".format(len(UNFOUND), ', '.join(UNFOUND)))

print_lcas(data/lca_result.json)
