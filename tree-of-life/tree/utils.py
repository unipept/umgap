import subprocess
from copy import deepcopy
from itertools import zip_longest

from tree.tree_of_life import Tree, Taxon


def print_tree_json(filename, tree, peptides, only_lcas=True):
    # Create a new tree
    reduced_tree = Tree()

    for peptide in peptides:
        for lineage in peptide.lineages:
            for taxon in lineage:
                if taxon and taxon != -1:
                    # Only add the taxon if it is not there yet
                    # Also don't add it when it's invalid and only_lcas is true
                    if not reduced_tree.taxons[taxon.taxon_id] and not (only_lcas and not taxon.valid_taxon):
                        reduced_tree.taxons[taxon.taxon_id] = Taxon(
                                taxon.taxon_id,
                                taxon.name,
                                taxon.rank,
                                taxon.get_parent(allow_no_rank=False, allow_invalid=True).taxon_id,
                                taxon.valid_taxon
                        )

                    # Break if we don't want to go deeper than the LCA
                    if only_lcas and taxon.taxon_id == peptide.lca:
                        break
            # Nifty little hack; also break out of this loop if we breaked
            # from the inner loop
            else:
                continue
            break

    reduced_tree.add_taxon_children()
    reduced_tree.add_taxon_parents()

    for peptide in peptides:
        if peptide.lca:
            reduced_tree.taxons[peptide.lca].self_count += 1

    reduced_tree.add_taxon_counts()
    reduced_tree.to_json(filename)


def compare_to_unipept(fastafile, pept2prot, tree, inputarray):

    unmatched = set()
    correct = 0

    # Get unipept LCAs
    unipept_lca_result = subprocess.Popen(
        "unipept pept2lca -i {} -s taxon_id -s peptide".format(fastafile),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    for lca in unipept_lca_result.stdout.readlines()[1:]:
        _, pept, lca = lca.decode('utf-8').strip().split(',')
        pept2prot[pept].unipept_lca = int(lca)


    for pept in inputarray:
        pept = pept2prot[pept]

        if pept.lca == pept.unipept_lca:
            correct += 1
        else:
            print("Unmatched for {}: {} ({}) vs {} ({})".format(pept.pept, tree.taxons[pept.lca].name, pept.lca, tree.taxons[pept.unipept_lca].name, pept.unipept_lca))
            unmatched.add(pept.pept)

    print("Unmatched: {}, {}".format(len(unmatched), ', '.join(unmatched)))
    print()
    print("Correct: {}, total: {}, accuracy: {}%".format(correct, len(inputarray), correct/len(inputarray)*100))
