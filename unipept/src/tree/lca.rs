
use std::ops::Add;
use std::collections::HashMap;

use agg::Aggregator;
use taxon;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;

pub struct LCACalculator {
    root: TaxonId,
    taxons: Vec<Option<Taxon>>,
    snapping: Vec<Option<TaxonId>>,
}

impl LCACalculator {
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self {
        let tree     = taxon::TaxonTree::new(&taxons);
        let by_id    = taxon::taxa_vector_by_id(taxons);
        let snapping = tree.snapping(&by_id, ranked_only);
        LCACalculator {
            root:     tree.root,
            taxons:   by_id,
            snapping: snapping,
        }
    }
}

impl Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        // Count the occurences
        let mut counts = HashMap::new();
        for taxon in taxons {
            *counts.entry(*taxon).or_insert(0) += 1;
        }

        // Build a subtree of the taxonomy
        let subtree = SubTree::new(self.root, &self.taxons, counts).unwrap().collapse(&Add::add);

        // TODO (for now, just an LCA implementation
        let lca = self.snapping[subtree.root].expect("LCA does not exist");
        self.taxons[lca].as_ref().expect("No LCA* found.")
    }
}
