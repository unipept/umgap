from functools import reduce
import os
import sys
import time

import numpy
import ctypes as ct

from lca_calculators.tree_of_life import get_tree, serialize_tree, RMQLIB_PATH
from lca_calculators.lca_calculator import LCA_Calculator
from lca_calculators.lca_calculator import DATA_DIR


class rmqinfo(ct.Structure):
    _fields_ = [
        ("alen", ct.c_int),
        ("array", ct.POINTER(ct.c_int)),
        ("sparse", ct.POINTER(ct.POINTER(ct.c_int))),
        ("block_min", ct.POINTER(ct.c_int)),
        ("labels", ct.POINTER(ct.c_int))
    ]


class Tree_LCA_Calculator(LCA_Calculator):

    def __init__(self):
        super().__init__()

        if os.path.isdir(DATA_DIR):
            starttime = time.time()
            print("Getting data from disk. Time: {}".format(starttime), file=sys.stderr)
            self.get_rmq_from_npy()
            self.get_tree_from_npy()
            print("Got data from disk: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
            print("---", file=sys.stderr)

        else:
            print("Warning: no seralized data found. Building tree from taxons.tsv and serializing the result. Use serialize=False to not store the result to disk.", file=sys.stderr)
            self.from_tree()
            self.store_data()
            self.get_tree_from_npy()


        self.librmq = ct.cdll.LoadLibrary(RMQLIB_PATH)
        self.preprocess_rmq()


    def store_data(self):
        starttime = time.time()
        print("Storing data to disk. Time: {}".format(starttime), file=sys.stderr)

        # Serialize the tree
        serialize_tree(self.tree)

        # Create the arrays
        if not os.path.isdir(DATA_DIR):
            os.makedirs(DATA_DIR)
        numpy.save(os.path.join(DATA_DIR, "euler_tour.npy"), self.euler_tour)
        numpy.save(os.path.join(DATA_DIR, "levels.npy"), self.levels)
        numpy.save(os.path.join(DATA_DIR, "first_occurences.npy"), self.first_occurences)

        print("Stored data to disk: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)


    def get_rmq_from_npy(self):
        """Gets the needed npy data from the preprocessed RMQ"""
        self.euler_tour = numpy.load(os.path.join(DATA_DIR, "euler_tour.npy"))
        self.levels = numpy.load(os.path.join(DATA_DIR, "levels.npy"))
        self.first_occurences = numpy.load(os.path.join(DATA_DIR, "first_occurences.npy"))


    def get_tree_from_npy(self):
        """Gets the needed npy data from the preprocessed tree"""
        self.names = numpy.load(os.path.join(DATA_DIR, "names.npy"))
        self.ranks = numpy.load(os.path.join(DATA_DIR, "ranks.npy"))
        self.valid_taxon_ids = numpy.load(os.path.join(DATA_DIR, "valid_taxon_ids.npy"))
        self.valid_ranked_taxon_ids = numpy.load(os.path.join(DATA_DIR, "valid_ranked_taxon_ids.npy"))


    def from_tree(self):
        self.tree = get_tree()

        starttime = time.time()
        print("Preprocessing LCA arrays. Time: {}".format(starttime), file=sys.stderr)

        self.euler_tour = numpy.array([None] * (2*self.tree.real_size - 1))
        self.levels = numpy.array([None] * (2*self.tree.real_size - 1))
        self.first_occurences = numpy.array([None] * self.tree.size)

        self.dfs_run(self.tree.taxons[1], 0, 0)

        print("Preprocessed LCA arrays: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)


    def preprocess_rmq(self):
        # Preprocess the RMQ
        starttime = time.time()
        print("Preprocessing RMQ. Time: {}".format(starttime), file=sys.stderr)

        self.librmq.rm_query_preprocess.restype = ct.POINTER(rmqinfo)
        self.rmqinfo = self.librmq.rm_query_preprocess(ct.c_void_p(self.levels.ctypes.data), len(self.levels))

        print("Preprocessed RMQ. Time: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
        print("---", file=sys.stderr)


    def free_rmqinfo(self):
        self.librmq.rm_free(self.rmqinfo)


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
        taxon_ids = (self.map_to_valid_taxon_id(taxon_id, allow_no_rank) for taxon_id in taxon_ids)

        #print("DEBUG:", file=sys.stderr)
        #print([self.tree.taxons[taxon].name for taxon in taxons], file=sys.stderr)

        # Put the highest join level as high as possible
        level_of_highest_join = sys.maxsize

        try:
            first = next(taxon_ids)
        except StopIteration:
            return None

        level_of_highest_join, taxon_id = reduce(self._calc_lca_pair, taxon_ids, (level_of_highest_join, first))

        return self.map_to_valid_taxon_id(taxon_id, allow_no_rank=allow_no_rank)


    def map_to_valid_taxon_id(self, taxon_id, allow_no_rank=False):
        if allow_no_rank:
            return self.valid_taxon_ids[taxon_id]
        else:
            return self.valid_ranked_taxon_ids[taxon_id]


    def get_rmq(self, start, end):
        return self.librmq.rm_query(self.rmqinfo, start, end)


    def cleanup(self):
        self.free_rmqinfo()
