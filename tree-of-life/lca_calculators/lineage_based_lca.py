from itertools import zip_longest
import time

from lca_calculators.lca_calculator import LCA_Calculator
from lca_calculators.tree_of_life import CLASSES


class Lineage_LCA_Calculator(LCA_Calculator):
    def __init__(self):
        super().__init__()

        self.genus_check = False


    def calc_lca(self, prots):
        """Given a list of protein ids, calculate the LCA"""
        lineages = []
        lca = None

        # Get the lineage for each prot
        for prot in prots:
            taxon = self.tree.taxons[int(prot)]
            lineage = taxon.get_lineage(allow_no_rank=False, allow_invalid=True)
            lineages.append(lineage)

        # If we have a result, get the LCA
        if lineages:
            lca = self.calc_lca_by_lineages(lineages, allow_invalid=False).taxon_id

        return lca


    def calc_lca_by_lineages(self, lineages, allow_invalid=True):
        """Does the actual LCA calculation"""

        # Use -1 as fillvalue here, we'll filter it out later
        for i, taxons in enumerate(zip_longest(*lineages, fillvalue=-1)):

            # Remove the filling and the invalid taxons if wanted
            if not allow_invalid:
                taxons = [t for t in taxons if t == -1 or t is None or t.valid_taxon]
            taxon_set = set(taxons) - set([-1])

            if self.genus_check:
                if i < len(CLASSES) and (CLASSES[i] == 'genus' or CLASSES[i] == 'species'):
                    taxon_set = taxon_set - set([None])

            if len(taxon_set) == 1:
                val = taxon_set.pop()
                if val:
                    lca = val
            elif len(taxon_set) > 1:
                return lca

        return lca


    def cleanup(self):
        pass
