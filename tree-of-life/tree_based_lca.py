import time

import vendor.rmq as rmq

class Tree_LCA_Calculator():

    def __init__(self, tree):
        self.euler_tour =[None] * (2*tree.real_size - 1)
        self.levels = [None] * (2*tree.real_size - 1)
        self.first_occurences = [None] * tree.size

        print("Preprocessing LCA arrays. Time: {}".format(time.time()))
        self.dfs_run(tree.taxons[1], 0, 0)
        print("Preprocessed LCA arrays: {}".format(time.time()))
        print("---")

        print("Preprocessing RMQ. Time: {}".format(time.time()))
        self.rmqinfo = rmq.rm_query_preprocess(self.levels, len(self.levels))
        print("Preprocessed RMQ. Time: {}".format(time.time()))
        print("---")


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


    def calc_lca(self, prots):
        if not prots:
            return None

        lca = prots[0]

        for prot in prots[1:]:
            lca_index = self.first_occurences[lca]
            prot_index = self.first_occurences[prot]
            #print("DEBUG: {}, {}:".format(lca_index, prot_index))
            rmq_index = self.get_rmq(lca_index, prot_index)

            # If one boundary is the result; take the other
            if rmq_index == lca_index:
                lca = self.euler_tour[prot_index]
            elif rmq_index == prot_index:
                lca = self.euler_tour[lca_index]
            else:
                lca = self.euler_tour[rmq_index]

        return lca


    def get_rmq(self, start, end):
        return rmq.rm_query(self.rmqinfo, start, end)
