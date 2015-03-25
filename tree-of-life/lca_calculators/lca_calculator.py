import time
import sys

from lca_calculators.tree_of_life import get_tree


class LCA_Calculator:
    """Provides methods to perform LCA calculations"""
    def __init__(self):
        print("Building tree. Time: {}".format(time.time()), file=sys.stderr)
        self.tree = get_tree()
        print("Tree built. Time: {}".format(time.time()), file=sys.stderr)
        print("---", file=sys.stderr)
