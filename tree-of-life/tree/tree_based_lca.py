import sys
import time

from tree.lca_calculator import LCA_Calculator
import tree.vendor.rmq as rmq


class Tree_LCA_Calculator(LCA_Calculator):

    def __init__(self):
        super().__init__()

        print("Preprocessing LCA arrays. Time: {}".format(time.time()), file=sys.stderr)
        self.euler_tour =[None] * (2*self.tree.real_size - 1)
        self.levels = [None] * (2*self.tree.real_size - 1)
        self.first_occurences = [None] * self.tree.size

        self.dfs_run(self.tree.taxons[1], 0, 0)

        print()
        print("Preprocessed LCA arrays: {}".format(time.time()), file=sys.stderr)
        print("---", file=sys.stderr)

        # Preprocess the RMQ
        print("Preprocessing RMQ. Time: {}".format(time.time()), file=sys.stderr)
        self.rmqinfo = self.preprocess_rmq()
        print("Preprocessed RMQ. Time: {}".format(time.time()), file=sys.stderr)
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


    def calc_lca(self, taxons):
        """Given a list of taxon ids, calculate the LCA"""

        if not taxons:
            return None

        # Map sort to their first valid parent
        taxons = [self.tree.taxons[taxon].get_parent(allow_no_rank=False, allow_invalid=False).taxon_id if not self.tree.taxons[taxon].valid_taxon else taxon for taxon in taxons]
        #print([self.tree.taxons[taxon].name for taxon in taxons], file=sys.stderr)

        # Don't do duplicate calculations, set the list first
        taxons = list(set(taxons))

        # Put the highest join level as high as possible
        level_of_highest_join = sys.maxsize

        # Iterate over the taxons, joining two taxons at a time
        while len(taxons) > 1:
            # Get two taxons
            first = taxons.pop()
            second = taxons.pop()

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

            taxons.append(self.euler_tour[lca_index])

            #print("Found LCA: {}".format(self.euler_tour[lca_index]), file=sys.stderr)
            #print("Level of highest join: {}".format(level_of_highest_join), file=sys.stderr)
            #print(file=sys.stderr)

        result = self.tree.taxons[taxons.pop()]

        if result.valid_taxon and result.rank != "no rank":
            return result.taxon_id
        else:
            return result.get_parent(allow_no_rank=False, allow_invalid=False).taxon_id


    def get_rmq(self, start, end):
        return rmq.rm_query(self.rmqinfo, start, end)


    def cleanup(self):
        self.free_rmqinfo()
