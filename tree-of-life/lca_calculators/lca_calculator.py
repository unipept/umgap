import time
import sys

from lca_calculators.tree_of_life import get_tree


class LCA_Calculator:
    """Provides methods to perform LCA calculations"""
    def __init__(self):
        starttime = time.time()
        print("Building tree. Time: {}".format(starttime), file=sys.stderr)
        self.tree = get_tree()
        print("Tree built. Time: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)
