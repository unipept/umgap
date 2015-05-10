# usage: cat peptides.fst | python commonpeptfilter.py pept2lca.fst

import sys


ranks_to_filter = set(["no rank", "superkingdom"])

lines = (line.strip().split(',') for line in open(sys.argv[1]))
peptides_to_filter = set(line[1] for line in lines if line[4] in ranks_to_filter)

for line in sys.stdin:
    if line.strip() not in peptides_to_filter:
        print(line.strip())
