"""Objects to hold a tree"""

import os
from json import JSONEncoder
import pickle

CLASSES = ['no rank', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus', 'species group', 'species subgroup', 'species', 'subspecies', 'varietas', 'forma']


class Tree():
    """It's all about the tree"""

    def __init__(self):
        self.size = 1500560 + 1
        self.taxons = [None] * self.size


    def from_taxons(self):
        """Read the tree from the taxons file"""
        self.read_taxons()
        self.add_taxon_children()
        self.add_taxon_parents()
        self.add_taxon_counts()

    def to_json(self, filename):
        """Writes the tree to JSON"""
        with open(filename, "wb") as f:
            self.taxons[1].to_json(f)


    def read_taxons(self):
        """Auxiliry method used to create a list of taxons"""

        def read_taxons_file():
            """Reads a taxon file and returns the tsv values as a splitted string"""
            with open('data/taxons.tsv') as taxon_file:
                # Skip the first two lines (`header` and `0 name 0 1`)
                next(taxon_file)
                next(taxon_file)

                for line in taxon_file:
                    yield line.split("\t")

        for line in read_taxons_file():
            self.taxons[int(line[0])] = Taxon(
                int(line[0]),
                line[1],
                line[2],
                int(line[3]),
                int(line[4])
            )


    def add_taxon_children(self):
        """Adds the children of each taxon"""
        for taxon in self.taxons:
            if taxon and taxon.taxon_id != taxon.parent_id:
                self.taxons[taxon.parent_id].children.add(taxon)

    def add_taxon_parents(self):
        """Adds parents to the taxons"""
        for taxon in self.taxons:
            if taxon:
                taxon.parent = self.taxons[taxon.parent_id]

    def add_taxon_counts(self):
        """Adds counts to the taxons"""
        self.taxons[1].add_counts()


class Taxon(JSONEncoder):
    """A Taxon objet"""

    def __init__(self, taxon_id, name, rank, parent_id, valid_taxon):
        super().__init__()

        self.taxon_id = taxon_id
        self.name = name
        self.rank = rank
        self.parent_id = parent_id
        self.valid_taxon = valid_taxon
        self.children = set()
        self.parent = None
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
"id":{self.taxon_id},
"name":"{self.name}",
"children":[
""".format(self=self).encode('utf-8'))

        for i, child in enumerate(self.children):
            child.to_json(f)
            if i != len(self.children) - 1: # JSON, y u no support trailing commas?
                f.write(",".encode('utf-8'))

        f.write("""
],
"data":{{"count":{self.count},"self_count":{self.self_count},"valid_taxon":{self.valid_taxon},"rank":"{self.rank}"}}}}
""".format(self=self).encode('utf-8'))


    def get_parent(self, allow_no_rank=True, allow_invalid=True):
        # Root case: return self
        if self.taxon_id == self.parent_id:
            return self

        if allow_no_rank or allow_invalid:
            return self.parent

        if not allow_no_rank and self.parent.rank != "no rank":
            return self.parent

        if not allow_invalid and not self.parent.valid_taxon:
            return self.parent

        return self.parent.get_parent(allow_no_rank, allow_invalid)


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
        if not allow_invalid and not self.valid_taxon:
            betweens = 0

        # Account for the no ranks, if wanted
        if not allow_no_rank and self.rank == "no rank":
            return parent_lineage

        return parent_lineage + between*betweens + [self]


def get_tree():
    """Creates or loads a tree object"""

    tree = Tree()
    tree.from_taxons()
    return tree
