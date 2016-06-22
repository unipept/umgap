
use std::ops::Add;

use agg::{Aggregator, count};
use taxon;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;

pub struct LCACalculator {
    pub root: TaxonId,
    pub taxons: Vec<Option<Taxon>>,
    pub snapping: Vec<Option<TaxonId>>,
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
        let counts  = count(taxons);
        let subtree = SubTree::new(self.root, &self.taxons, counts).unwrap().collapse(&Add::add);
        let lca     = self.snapping[subtree.root].expect("LCA does not exist");
        self.taxons[lca].as_ref().expect("No LCA* found.")
    }
}
