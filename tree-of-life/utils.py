import subprocess
from itertools import zip_longest

def print_tree_json(filename, tree, lcas):
    lineages = [tree.taxons[prot].get_lineage(no_rank=False, invalid=False) for prot in lcas if prot]

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
            taxon.self_count = lcas.count(taxon.taxon_id)

    # Percolate the counts up the tree
    lineages[0][0].add_counts()

    # Print the result
    with open(filename, "wb") as f:
        lineages[0][0].to_json(f)


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
