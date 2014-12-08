"""Calculates the LCA for a given protein"""

import sys
import subprocess
from itertools import zip_longest

from tree_of_life import get_tree, CLASSES

TREE = get_tree()
LCAS = []
UNFOUND = set()
UNMATCHED = set()
CORRECT = 0

def print_tree_json(filename):
    lineages = [TREE.taxons[prot].get_lineage(no_rank=False, invalid=False) for prot in LCAS if prot]

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
            taxon.self_count = LCAS.count(taxon.taxon_id)

    # Percolate the counts up the tree
    lineages[0][0].add_counts()

    # Print the result
    with open(filename, "wb") as f:
        lineages[0][0].to_json(f)


class Peptide:
    def __init__(self, pept):
        self.pept = pept
        self.prots = []
        self.lca = None
        self.unipept_lca = None
        self.lineages = []


    def addprot(self, prot):
        self.prots.append(prot)


    def handle_lca(self):
        """Gets the lineages from the tree and performs the LCA calculation"""

        # Get the lineage for each prot
        for prot in self.prots:
            taxon = TREE.taxons[int(prot)]
            lineage = taxon.get_lineage(no_rank=False, invalid=False)
            self.lineages.append(lineage)

        # If we have a result, get the LCA, otherwise, add it to the unfound ones
        if self.lineages:
            self.lca = self.get_lca().taxon_id
        else:
            UNFOUND.add(self.pept)


    def get_lca(self):
        """Does the actual LCA calculation"""

        # Use -1 as fillvalue here, we'll filter it out later
        for i, taxons in enumerate(zip_longest(*self.lineages, fillvalue=-1)):
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


fastafile = "/tmp/fasta.tmp"
# Dict containing the input sequence
inputarray = []
# Dict containing the peptides to their objects
pept2prot = dict()

# Create a temp fastafile for easy Unipept querying
with open(fastafile, "wb") as f:
    for i, line in enumerate(sys.stdin):
        f.write(">|{}\n".format(i).encode('utf-8'))
        f.write(line.encode('utf-8'))

        line = line.strip()

        inputarray.append(line)

        if line not in pept2prot:
            pept2prot[line] = Peptide(line.strip())


# Get all the proteins
prot_result = subprocess.Popen(
    "unipept pept2prot -i {} -s taxon_id -s peptide".format(fastafile),
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

for prot in prot_result.stdout.readlines()[1:]:
    _, pept, prot = prot.decode('utf-8').strip().split(',')
    pept2prot[pept].addprot(prot)


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

# Get all the LCAs
for pept in pept2prot.values():
    pept.handle_lca()


for pept in inputarray:
    pept = pept2prot[pept]

    if pept.lca == pept.unipept_lca:
        CORRECT += 1
        LCAS.append(pept.lca)
    else:
        UNMATCHED.add(pept.pept)

print("Unmatched: {}, {}".format(len(UNMATCHED), ', '.join(UNMATCHED)))
print("Unfound: {}, {}".format(len(UNFOUND), ', '.join(UNFOUND)))
print()
print("Correct: {}, total: {}, accuracy: {}%".format(CORRECT, len(inputarray), CORRECT/len(inputarray)*100))

print_tree_json("data/lca_result.json")
