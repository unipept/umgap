from tree_of_life import get_tree
t = get_tree()
t.from_taxons()
t.to_json("hierarchical-taxons.json")
