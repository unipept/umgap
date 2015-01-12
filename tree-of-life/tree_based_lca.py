import time

import vendor.rmq as rmq

class Real_LCA_Calculator():

    def __init__(self, tree):
        self.euler_tour =[None] * (2*tree.real_size - 1)
        self.levels = [None] * (2*tree.real_size - 1)
        self.first_occurences = [None] * tree.size

        print("Preprocessing LCA arrays. Time: {}".format(time.time()))
        self.dfs_run(tree.taxons[1], 0, 0)
        print("Preprocessed LCA arrays: {}".format(time.time()))

        print("Preprocessing RMQ. Time: {}".format(time.time()))
        self.rmqinfo = rmq.rm_query_preprocess(self.levels, len(self.levels))
        print("Preprocessed RMQ. Time: {}".format(time.time()))

        print("RMQ Test 0:10000. Time: {}".format(time.time()))
        self.get_rmq(10000, 15000)
        print("RMQ Tested 0:10000. Time: {}".format(time.time()))


    def dfs_run(self, taxon, iteration, level):
        """Method to fill the RMQ arrays, we walk DFS trough the tree"""

        self.euler_tour[iteration] = taxon.taxon_id
        self.levels[iteration] = level
        if not self.first_occurences[taxon.taxon_id]:
            self.first_occurences[taxon.taxon_id] = iteration

        for child in taxon.children:
            iteration = self.dfs_run(child, iteration + 1, level + 1)

            self.euler_tour[iteration] = taxon.taxon_id
            self.levels[iteration] = level

        return iteration + 1

    def get_rmq(self, start, end):
        print(rmq.rm_query(self.rmqinfo, start, end))
