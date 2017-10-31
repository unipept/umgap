"""Calculates the LCA from stdin given a list of proteins"""

import sys
import subprocess
import time
from collections import OrderedDict

from lca_calculators.tree_based_lca import Tree_LCA_Calculator
from lca_calculators.lineage_based_lca import Lineage_LCA_Calculator
from lca_calculators.utils import print_tree_json, compare_to_unipept


fastafile = "/tmp/fasta.tmp"
# Dict containing the input sequence
inputarray = []
# Dict containing the peptides to their objects
pept2prot = OrderedDict() # Uncomment this to enable diffing between resultsets
#pept2prot = dict()

lcas = []
unfound = set()

class Peptide:
    def __init__(self, pept):
        self.pept = pept
        self.prots = []
        self.lca = None
        self.unipept_lca = None


    def add_prot(self, prot):
        self.prots.append(prot)


# Create a temp fastafile for easy Unipept querying
print("Create fasta. Time: {}".format(time.time()), file=sys.stderr)
with open(fastafile, "wb") as f:
    for i, line in enumerate(sys.stdin):
        f.write(">|{}\n".format(i).encode('utf-8'))
        f.write(line.encode('utf-8'))

        line = line.strip()

        inputarray.append(line)

        if line not in pept2prot:
            pept2prot[line] = Peptide(line.strip())

print("Created fasta. Time: {}".format(time.time()), file=sys.stderr)
print("---", file=sys.stderr)


# Get all the proteins
print("Get proteins. Time: {}".format(time.time()), file=sys.stderr)
prot_result = subprocess.Popen(
    "unipept pept2prot -i {} -s taxon_id -s peptide".format(fastafile),
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

# Map every peptide to a list of proteins
for prot in prot_result.stdout.readlines()[1:]:
    _, pept, prot = prot.decode('utf-8').strip().split(',')

    pept2prot[pept].add_prot(int(prot))

print("Got proteins. Time: {}".format(time.time()), file=sys.stderr)
print("---", file=sys.stderr)

# Calculate all the LCAs
calculator = Tree_LCA_Calculator()
#calculator = Lineage_LCA_Calculator()
print("Get Tree LCAs. Time: {}".format(time.time()), file=sys.stderr)
for pept in pept2prot.values():
    print("Calc for {}".format(pept.pept))
    pept.lca = calculator.calc_lca(pept.prots)

    if pept.lca:
        print("LCA for {} is {} ({})".format(pept.pept, calculator.tree.taxons[pept.lca].name, pept.lca))
        lcas.append(pept.lca)
    else:
        print("LCA for {} not found".format(pept.pept))
        unfound.add(pept)

calculator.cleanup()
print("Got Tree LCAs. Time: {}".format(time.time()), file=sys.stderr)


#if UNFOUND:
#    print("Unfound: {}, {}".format(len(UNFOUND), ', '.join(UNFOUND)))

compare_to_unipept(fastafile, pept2prot, calculator.tree, inputarray)
#print_tree_json("visualisation-data/lca_with_overall_check.json", TREE, pept2prot.values(), only_lcas=True)
