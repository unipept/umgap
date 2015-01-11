"""Calculates the LCA for a given protein"""

import sys
import subprocess
from itertools import zip_longest
from collections import OrderedDict

from tree_of_life import get_tree, CLASSES
from utils import print_tree_json, compare_to_unipept

GENUS_CHECK = True
DEBUG = False
TREE = get_tree()
LCAS = []
UNFOUND = set()


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
            lineage = taxon.get_lineage(allow_no_rank=False, allow_invalid=True)
            self.lineages.append(lineage)

        if DEBUG:
            for lineage in self.lineages:
                print([(taxon.taxon_id, taxon.valid_taxon) if taxon and taxon != -1 else -1 if taxon == -1 else 0 for taxon in lineage])
            print()

        # If we have a result, get the LCA, otherwise, add it to the unfound ones
        if self.lineages:
            self.lca = self.get_lca(allow_invalid=False).taxon_id
            LCAS.append(self.lca)
            print("LCA for {} is {} ({})".format(self.pept, TREE.taxons[self.lca].name, self.lca))
        else:
            UNFOUND.add(self.pept)
            print("LCA for {} not found".format(self.pept))



    def get_lca(self, allow_invalid=True):
        """Does the actual LCA calculation"""

        # Use -1 as fillvalue here, we'll filter it out later
        for i, taxons in enumerate(zip_longest(*self.lineages, fillvalue=-1)):

            if DEBUG:
                print("Before: {}".format([(taxon.taxon_id, taxon.valid_taxon) if taxon and taxon != -1 else -1 if taxon == -1 else 0 or taxon in taxons]))

            # Remove the filling and the invalid taxons if wanted
            if not allow_invalid:
                taxons = [t for t in taxons if t == -1 or t is None or t.valid_taxon]
            taxon_set = set(taxons) - set([-1])

            if DEBUG:
                print("After: {}".format([(taxon.taxon_id, taxon.valid_taxon) if taxon and taxon != -1 else -1 if taxon == -1 else 0 or taxon in taxons]))

            if GENUS_CHECK:
                if i < len(CLASSES) and (CLASSES[i] == 'genus' or CLASSES[i] == 'species'):
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
#pept2prot = OrderedDict() # Uncomment this to enable diffing between resultsets
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

# Get all the LCAs
for pept in pept2prot.values():
    pept.handle_lca()

#if UNFOUND:
#    print("Unfound: {}, {}".format(len(UNFOUND), ', '.join(UNFOUND)))

#compare_to_unipept(fastafile, pept2prot, TREE, inputarray)
#print_tree_json("visualisation-data/lca_with_overall_check.json", TREE, pept2prot.values(), only_lcas=True)
