"""Turns a list of taxons into a visualisation json
Inputs can be:
```
    taxon1
    taxon2
    taxon1
    taxon1
    taxon2
    ...
```

or

```
    taxon1,count1
    taxon2,count2
```"""

import sys

from lca_calculators.tree_of_life import Tree, Taxon, get_tree

ALLOW_NO_RANK=False
ALLOW_INVALID=False
FILE=sys.argv[1]

# Get the tree
tree = get_tree()

taxons = dict()

for line in sys.stdin:

    line = line.strip()

    if ',' in line:
        taxon_id, count = line.split(',')
    else:
        taxon_id = line
        count = 1

    taxon_id = int(taxon_id)
    count = int(count)

    if taxon_id not in taxons:
        taxons[taxon_id] = count
    else:
        taxons[taxon_id] += count


reduced_tree = Tree()

for lineage in [tree.taxons[taxon].get_lineage(allow_no_rank=ALLOW_NO_RANK, allow_invalid=ALLOW_INVALID) for taxon in taxons]:
    for taxon in lineage:
        if taxon and taxon != -1:
            # Only add the taxon if it is not there yet
            # Also don't add it when it's invalid and only_lcas is true
            if not reduced_tree.taxons[taxon.taxon_id]:
                reduced_tree.taxons[taxon.taxon_id] = Taxon(
                        taxon.taxon_id,
                        taxon.name,
                        taxon.rank,
                        taxon.get_parent(allow_no_rank=False, allow_invalid=True).taxon_id,
                        taxon.valid_taxon
                )


reduced_tree.add_children_to_taxons()
reduced_tree.add_parents_to_taxons()

for taxon in taxons:
    reduced_tree.taxons[taxon].self_count = taxons[taxon]

reduced_tree.add_counts_to_taxons()
reduced_tree.to_json(FILE)
