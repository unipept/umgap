# usage: cat peptides.fst | python commonpeptfilter.py pept2lca.fst

import sys


ranks_to_filter = set(["no rank", "superkingdom"])

lines = ((line.strip().split(',')) for line in open(sys.argv[1]))

# {pept: [pept2lca]}
peptides_to_filter = { line[1]: line for line in lines if line[4] in ranks_to_filter }

with open(sys.argv[2], 'w') as f:
    for pept in sys.stdin:
        pept = pept.strip()

        if pept not in peptides_to_filter:
            print(pept)
        else:
            line = peptides_to_filter[pept]
            f.write("{},{},,{}\n".format(line[0], line[1], line[2]))
