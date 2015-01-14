import time

import tree_of_life


class LCA_Calculator:
    """Provides methods to perform LCA calculations"""
    def __init__(self):
        print("Building tree. Time: {}".format(time.time()))
        self.tree = tree_of_life.get_tree()
        print("Tree built. Time: {}".format(time.time()))
        print("---")
