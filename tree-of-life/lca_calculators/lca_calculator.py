import time
import sys

from lca_calculators.tree_of_life import get_tree


class LCA_Calculator:
    """Provides methods to perform LCA calculations"""
    def __init__(self):
        self.tree = get_tree()
