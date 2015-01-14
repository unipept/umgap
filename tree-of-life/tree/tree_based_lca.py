import sys
import time

from lca_calculator import LCA_Calculator
import vendor.rmq as rmq


class Tree_LCA_Calculator(LCA_Calculator):

    def __init__(self):
        super().__init__()

        self.euler_tour =[None] * (2*self.tree.real_size - 1)
        self.levels = [None] * (2*self.tree.real_size - 1)
        self.first_occurences = [None] * self.tree.size

        print("Preprocessing LCA arrays. Time: {}".format(time.time()))
        self.dfs_run(self.tree.taxons[1], 0, 0)
        print("Preprocessed LCA arrays: {}".format(time.time()))
        print("---")

        # Preprocess the RMQ
        print("Preprocessing RMQ. Time: {}".format(time.time()))
        self.preprocess_rmq()
        print("Preprocessed RMQ. Time: {}".format(time.time()))
        print("---")


    def preprocess_rmq(self):
        self.rmqinfo = rmq.rm_query_preprocess(self.levels, len(self.levels))


    def free_rmqinfo(self):
        rmq.rm_free(self.free_rmqinfo)


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
        """Given a list of protein ids, calculate the LCA"""

        if not prots:
            return None

        # Map prots to their first valid parent
        prots = [TREE.taxons[prot].get_parent(allow_invalid=False).taxon_id for prot in self.prots]

        lca = prots[0]

        for prot in prots[1:]:
            lca_index = self.first_occurences[lca]
            prot_index = self.first_occurences[prot]
            #print("DEBUG: {}, {}, {}, {}:".format(lca_index, prot_index, lca_index < len(self.levels), prot_index < len(self.levels)))
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


    def cleanup(self):
        self.free_rmqinfo()
