from lca_calculators.tree_based_lca import Tree_LCA_Calculator
calc = Tree_LCA_Calculator()

taxons = [7955, 41, 54]
print(calc.calc_lca(taxons))
