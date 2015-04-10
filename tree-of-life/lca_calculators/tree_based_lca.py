from functools import reduce
import sys
import time

from lca_calculators.lca_calculator import LCA_Calculator
import lca_calculators.vendor.rmq as rmq


class Tree_LCA_Calculator(LCA_Calculator):

    def __init__(self):
        super().__init__()

        starttime = time.time()
        print("Preprocessing LCA arrays. Time: {}".format(starttime), file=sys.stderr)
        self.euler_tour =[None] * (2*self.tree.real_size - 1)
        self.levels = [None] * (2*self.tree.real_size - 1)
        self.first_occurences = [None] * self.tree.size

        self.dfs_run(self.tree.taxons[1], 0, 0)

        print()
        print("Preprocessed LCA arrays: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)

        # Preprocess the RMQ
        starttime = time.time()
        print("Preprocessing RMQ. Time: {}".format(starttime), file=sys.stderr)
        self.rmqinfo = self.preprocess_rmq()
        print("Preprocessed RMQ. Time: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)


    def preprocess_rmq(self):
        return rmq.rm_query_preprocess(self.levels, len(self.levels))


    def free_rmqinfo(self):
        rmq.rm_free(self.rmqinfo)


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


    def _calc_lca_pair(self, acc, second):
        level_of_highest_join, first = acc

        # Get their occurence in the levels/euler_tour array
        first_index = self.first_occurences[first]
        second_index = self.first_occurences[second]

        #print("DEBUG:", file=sys.stderr)
        #print("LCA BY {} AND {} ".format(
            #self.tree.taxons[self.euler_tour[first_index]].taxon_id,
            #self.tree.taxons[self.euler_tour[second_index]].taxon_id
        #), file=sys.stderr)

        # Calculate the index of their LCA
        rmq_index = self.get_rmq(first_index, second_index)

        #print("LCA BY IS {} AT LEVEL {} ".format(
            #self.tree.taxons[self.euler_tour[rmq_index]].taxon_id,
            #self.levels[rmq_index]
        #), file=sys.stderr)

        # We've found a split, take it and update the highest_join level
        if rmq_index != first_index and rmq_index != second_index:
            lca_index = rmq_index
            level_of_highest_join = self.levels[rmq_index]
        else:
            # We're on a single lineage here, take the bottom one
            if rmq_index == first_index:
                lca_index = second_index
            elif rmq_index == second_index:
                lca_index = first_index

            # Don't go lower than the max join though
            if self.levels[lca_index] > level_of_highest_join:
                lca_index = rmq_index


        #print("Found LCA: {}".format(self.euler_tour[lca_index]), file=sys.stderr)
        #print("Level of highest join: {}".format(level_of_highest_join), file=sys.stderr)
        #print(file=sys.stderr)

        return level_of_highest_join, self.euler_tour[lca_index]


    def calc_lca(self, taxon_ids, allow_no_rank=False):
        """Given a list of taxon ids, calculate the LCA"""

        # Map taxons their valid counterpart
        taxon_ids = (self.tree.taxons[taxon_id].map_to_valid_taxon_id(allow_no_rank=allow_no_rank)
                  for taxon_id in taxon_ids)

        #print("DEBUG:", file=sys.stderr)
        #print([self.tree.taxons[taxon].name for taxon in taxons], file=sys.stderr)

        # Put the highest join level as high as possible
        level_of_highest_join = sys.maxsize

        try:
            first = next(taxon_ids)
        except StopIteration:
            return None

        level_of_higest_join, taxon_id = reduce(self._calc_lca_pair, taxon_ids, (level_of_highest_join, first))
        taxon = self.tree.taxons[taxon_id]

        return taxon.map_to_valid_taxon_id(allow_no_rank=allow_no_rank)


    def get_rmq(self, start, end):
        return rmq.rm_query(self.rmqinfo, start, end)


    def cleanup(self):
        self.free_rmqinfo()
