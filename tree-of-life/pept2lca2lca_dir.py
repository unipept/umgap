"""Expects an input sorted by fasta ID"""

import sys
import glob
import os


from tree.tree_based_lca import Tree_LCA_Calculator

calculator = Tree_LCA_Calculator()

os.chdir(argv[1])

print("fasta_header,taxon_id,taxon_name,taxon_rank")
for file in glob.glob("*.fst"):
    with open(file, 'r') as f:

        # Skip the header
        next(sys.stdin)

        ids = []

        for line in f.readlines():
            splitted_line = line.strip().split(',')

            fasta_id = splitted_line[0][2:]
            ids.append(splited_line[3])

        lca = calculator.calc_lca(int(ids))

        taxon = calculator.tree.taxons[lca]
        print(">|{},{},{},{}".format(fasta_id, taxon.taxon_id, taxon.name, taxon.rank))
