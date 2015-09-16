"""Objects to hold a tree"""

from json import JSONEncoder
import os
import sys
import time

import numpy

CLASSES = ['no rank', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus', 'species group', 'species subgroup', 'species', 'subspecies', 'varietas', 'forma']

RMQLIB_PATH = os.path.join(os.path.dirname(__file__), 'vendor/librmq.o')
DEFAULT_DATA_DIR = os.path.join(os.path.expanduser('~'), ".rmqnpydata")


class Tree():
    """It's all about the tree"""

    def __init__(self):
        with open(os.path.join(os.path.dirname(__file__), 'data/taxons.tsv')) as taxon_file:
            for i, l in enumerate(taxon_file):
                pass

            self.size = int(l.split()[0])+1
            self.real_size = i+1
            self.taxons = [None] * self.size


    def from_taxons(self):
        """Read the tree from the taxons file"""
        self.read_taxons()

        # We cant put these two in the loop as all the parents and children
        # need to have been added first
        self.add_children_to_taxons()
        self.add_parents_to_taxons()

        for taxon in self.taxons:
            self.add_valid_parent_to_taxon(taxon)
            self.add_valid_ranked_parent_to_taxon(taxon)
            self.add_valid_taxon_to_taxon(taxon)
            self.add_valid_ranked_taxon_to_taxon(taxon)

        self.add_counts_to_taxons()


    def to_json(self, filename):
        """Writes the tree to JSON"""
        with open(filename, "wb") as f:
            self.taxons[1].to_json(f)


    def read_taxons(self):
        """Auxiliary method used to create a list of taxons"""

        def read_taxons_file():
            """Reads a taxon file and returns the tsv values as a splitted string"""
            with open(os.path.join(os.path.dirname(__file__), 'data/taxons.tsv')) as taxon_file:
                for line in taxon_file:
                    yield line.rstrip().split("\t")


        for line in read_taxons_file():
            self.taxons[int(line[0])] = Taxon(
                int(line[0]),
                line[1],
                line[2],
                int(line[3]),
                '1' == format(ord(line[4][0]), 'b')
            )


    def add_children_to_taxons(self):
        """Adds the children of each taxon"""
        for taxon in self.taxons:
            if taxon and taxon.taxon_id != taxon.parent_id:
                self.taxons[taxon.parent_id].children.add(taxon)


    def add_parents_to_taxons(self):
        """Adds parents to the taxons"""
        for taxon in self.taxons:
            if taxon:
                taxon.parent = self.taxons[taxon.parent_id]


    def add_valid_parent_to_taxon(self, taxon):
        """Adds a valid_parent_id to the taxons"""
        if taxon:
            taxon.valid_parent_id = taxon.get_parent(allow_invalid=False).taxon_id


    def add_valid_ranked_parent_to_taxon(self, taxon):
        """Adds a ranked_valid_parent_id to the taxons"""
        if taxon:
            taxon.valid_ranked_parent_id = taxon.get_parent(allow_no_rank=False, allow_invalid=False).taxon_id


    def add_valid_taxon_to_taxon(self, taxon):
        """Adds a valid_parent_id to the taxons"""
        if taxon:
            taxon.valid_taxon_id = taxon.map_to_valid_taxon_id()


    def add_valid_ranked_taxon_to_taxon(self, taxon):
        """Adds a valid_parent_id to the taxons"""
        if taxon:
            taxon.valid_ranked_taxon_id = taxon.map_to_valid_taxon_id(allow_no_rank=False)



    def add_counts_to_taxons(self):
        """Adds counts to the taxons"""
        self.taxons[1].add_counts()


class Taxon(JSONEncoder):
    """A Taxon objet"""

    def __init__(self, taxon_id, name, rank, parent_id, valid_taxon):
        super().__init__()

        self.taxon_id = taxon_id
        self.name = name
        self.rank = rank

        self.parent = None
        self.children = set()

        self.parent_id = parent_id
        self.valid_parent_id = None
        self.valid_ranked_parent_id = None
        self.valid_taxon_id = None
        self.valid_ranked_taxon_id = None

        self.valid_taxon = valid_taxon

        self.count = 0
        self.self_count = 0

    def add_counts(self):
        children = self.self_count

        for child in self.children:
            children += child.add_counts()

        self.count = children

        return children

    def add_self_counts(self):
        self.self_count = 0 if len(self.children) else 1


    def to_json(self, f):
        """This must be the worst printer I've ever written."""
        f.write("""
{{
"id":{taxon_id},
"name":"{name}",
"children":[
""".format(
        taxon_id=self.taxon_id,
        name=self.name,
    ).encode('utf-8'))

        for i, child in enumerate(self.children):
            child.to_json(f)
            if i != len(self.children) - 1: # JSON, y u no support trailing commas?
                f.write(",".encode('utf-8'))

        f.write("""
],
"data":{{"count":{count},"self_count":{self_count},"valid_taxon":{valid_taxon},"rank":"{rank}"}}}}
""".format(
        count=self.count,
        self_count=self.self_count,
        valid_taxon=int(self.valid_taxon),
        rank=self.rank,
    ).encode('utf-8'))


    def get_parent(self, allow_no_rank=True, allow_invalid=True):
        # Root case: return self
        if self.taxon_id == self.parent_id:
            return self

        # A big if else thing, makes it easier to reason about it
        if allow_no_rank:
            if allow_invalid:
                # Allow no ranks and invalids: return all
                return self.parent
            else:
                # Allow no ranks but no invalids
                if self.parent.valid_taxon:
                    # If the parent is valid: return the parent
                    return self.parent
                else:
                    # If the parent is invalid, return its valid parent
                    return self.parent.get_parent(allow_no_rank, allow_invalid)
        else:
            if allow_invalid:
                # Dont allow no ranks, allow invalids
                if self.parent.rank != "no rank":
                    # If the parent is no no rank, return it
                    return self.parent
                else:
                    # If the parent is a no rank, return its valid parent
                    return self.parent.get_parent(allow_no_rank, allow_invalid)
            else:
                # Dont allow no ranks, don't allow invalids
                if self.parent.rank != "no rank" and self.parent.valid_taxon:
                    return self.parent
                else:
                    return self.parent.get_parent(allow_no_rank, allow_invalid)


    def is_ranked(self):
        return self.rank != "no rank"


    def map_to_valid_taxon_id(self, allow_no_rank=True):
        if not self.valid_taxon:
            if allow_no_rank:
                return self.valid_parent_id
            else:
                return self.valid_ranked_parent_id
        else:
            if allow_no_rank or self.is_ranked():
                return self.taxon_id
            else:
                return self.valid_ranked_parent_id


    def get_lineage(self, allow_no_rank=True, allow_invalid=True):
        """Gets the lineage of a taxon to the root"""

        # Root, base case
        if self.taxon_id == self.parent_id:
            return [self]

        # Calculate the distance between a child and its parent
        between = [None]
        betweens = 0
        parent = self.get_parent(allow_no_rank, allow_invalid)
        parent_lineage = parent.get_lineage(allow_no_rank, allow_invalid)

        parent_rank_id = CLASSES.index(parent.rank)
        self_rank_id = CLASSES.index(self.rank)
        betweens = self_rank_id - parent_rank_id - 1

        # Account for the invalids, if wanted
        if not self.valid_taxon:
            if allow_invalid:
                between = [-1]
            else:
                return parent_lineage

        # Account for the no ranks, if wanted
        if not allow_no_rank and self.rank == "no rank":
            return parent_lineage

        return parent_lineage + between*betweens + [self]


def get_tree():
    """Creates or loads a tree object"""

    starttime = time.time()
    print("Building tree. Time: {}".format(starttime), file=sys.stderr)

    tree = Tree()
    tree.from_taxons()

    print("Built tree. Time: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
    print("---", file=sys.stderr)

    return tree


def serialize_tree(tree=None, rmqdatadir=DEFAULT_DATA_DIR):
    if not tree:
        tree = get_tree()

    starttime = time.time()
    print("Serialize tree. Time: {}".format(starttime), file=sys.stderr)

    names = numpy.array([None] * tree.size)
    ranks = numpy.array([None] * tree.size)
    valids = numpy.array([None] * tree.size)
    valid_taxon_ids = numpy.array([None] * tree.size)
    valid_ranked_taxon_ids = numpy.array([None] * tree.size)

    for taxon in tree.taxons:
        if taxon:
            names[taxon.taxon_id]                   = taxon.name
            ranks[taxon.taxon_id]                   = taxon.rank
            valids[taxon.taxon_id]                  = taxon.valid_taxon
            valid_taxon_ids[taxon.taxon_id]         = taxon.valid_taxon_id
            valid_ranked_taxon_ids[taxon.taxon_id]  = taxon.valid_ranked_taxon_id

    if not os.path.isdir(rmqdatadir):
        os.makedirs(rmqdatadir)
    numpy.save(os.path.join(rmqdatadir, "names.npy"), names)
    numpy.save(os.path.join(rmqdatadir, "ranks.npy"), ranks)
    numpy.save(os.path.join(rmqdatadir, "valids.npy"), valids)
    numpy.save(os.path.join(rmqdatadir, "valid_taxon_ids.npy"), valid_taxon_ids)
    numpy.save(os.path.join(rmqdatadir, "valid_ranked_taxon_ids.npy"), valid_ranked_taxon_ids)

    print("Serialized tree. Time: {}, time elapsed: {}".format(time.time(), time.time()-starttime), file=sys.stderr)
    print("---", file=sys.stderr)
