"""Objects to hold a tree"""

import os
from json import JSONEncoder
import pickle

class Tree():
    """It's all about the tree"""

    def __init__(self):
        self.size = 1500560 + 1
        self.taxons = [None] * self.size

        self.read_taxons()
        self.add_taxon_children()
        self.add_taxon_parents()


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
        for taxon in self.taxons:
            if taxon:
                taxon.parent = self.taxons[taxon.parent_id]


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

    def to_json(self, f):
        """This must be the worst printer I've ever written."""
        f.write("""
{{
"id":{self.taxon_id},
"name":"{self.name}",
"children":[
""".format(self=self).encode('utf-8'))

        children = 1

        for i, child in enumerate(self.children):
            children += child.to_json(f)
            if i != len(self.children) - 1: # JSON, y u no support trailing commas?
                f.write(",".encode('utf-8'))

        f.write("""
],
"data":{{"count":{children},"self_count":{self_count},"valid_taxon":{self.valid_taxon},"rank":"{self.rank}"}}}}
""".format(children=children, self_count=len(self.children)+1, self=self).encode('utf-8'))

        return children

    def get_lineage(self, no_rank=True, invalid=True):
        """Gets the lineage of a taxon to the root"""

        # Root, return an empty list
        if self.taxon_id == self.parent_id:
            return [self.taxon_id]

        # Don't count no ranks if we don't want to
        if not no_rank and self.rank == "no_rank":
            return self.parent.get_lineage()

        # Don't count invalids if we don't want to
        if not invalid and not self.valid_taxon:
            return self.parent.get_lineage()

        return self.parent.get_lineage() + [self.taxon_id]


def get_tree():
    """Creates or loads a tree object"""

    # Return the tree if it exists
    # if os.path.exists('tree.p'):
    #     print("Loading tree from pickle")
    #     tree = pickle.load(open("tree.p", "rb"))
    # else: # Create it otherwise
    tree = Tree()
    #     pickle.dump(tree, open("tree.p", "wb"))

    return tree
