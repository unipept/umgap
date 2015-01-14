"""Calculates the LCA from stdin given a list of proteins"""

import sys
import subprocess
import time
from collections import OrderedDict

from tree_based_lca import Tree_LCA_Calculator
from lineage_based_lca import Lineage_LCA_Calculator
from utils import print_tree_json, compare_to_unipept


fastafile = "/tmp/fasta.tmp"
# Dict containing the input sequence
inputarray = []
# Dict containing the peptides to their objects
#pept2prot = OrderedDict() # Uncomment this to enable diffing between resultsets
pept2prot = dict()

lcas = []
unfound = set()


# Create a temp fastafile for easy Unipept querying
print("Create fasta. Time: {}".format(time.time()))
with open(fastafile, "wb") as f:
    for i, line in enumerate(sys.stdin):
        f.write(">|{}\n".format(i).encode('utf-8'))
        f.write(line.encode('utf-8'))

        line = line.strip()

        inputarray.append(line)
print("Created fasta. Time: {}".format(time.time()))
print("---")


# Get all the proteins
print("Get proteins. Time: {}".format(time.time()))
prot_result = subprocess.Popen(
    "unipept pept2prot -i {} -s taxon_id -s peptide".format(fastafile),
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT
)

# Map every peptide to a list of proteins
for prot in prot_result.stdout.readlines()[1:]:
    _, pept, prot = prot.decode('utf-8').strip().split(',')

    if pept in pept2prot:
        pept2prot[pept].append(prot)
    else:
        pept2prot[pept] = [prot]
print("Got proteins. Time: {}".format(time.time()))
print("---")

# Calculate all the LCAs
print("Get Tree LCAs. Time: {}".format(time.time()))
calculator = Lineage_LCA_Calculator()
for pept in pept2prot:
    lca = calculator.calc_lca(pept2prot[pept])

    if lca:
        print("LCA for {} is {} ({})".format(pept, calculator.tree.taxons[lca].name, lca))
        lcas.append(lca)
    else:
        print("LCA for {} not found".format(pept))
        unfound.add(pept)
calculator.cleanup()
print("Got Tree LCAs. Time: {}".format(time.time()))


#if UNFOUND:
#    print("Unfound: {}, {}".format(len(UNFOUND), ', '.join(UNFOUND)))

#compare_to_unipept(fastafile, pept2prot, TREE, inputarray)
#print_tree_json("visualisation-data/lca_with_overall_check.json", TREE, pept2prot.values(), only_lcas=True)
