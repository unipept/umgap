
use std::ops::Add;

use agg::{Aggregator, count};
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;

pub struct LCACalculator {
    pub root: TaxonId,
    pub parents: Vec<Option<TaxonId>>,
}

impl LCACalculator {
    pub fn new(root: TaxonId, taxonomy: &Vec<Option<Taxon>>) -> Self {
        LCACalculator {
            root:    root,
            parents: taxonomy.iter().map(|mt| mt.as_ref().map(|t| t.parent)).collect()
        }
    }
}

impl Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> TaxonId {
        let counts  = count(taxons);
        let subtree = SubTree::new(self.root, &self.parents, counts).unwrap().collapse(&Add::add);
        subtree.root
    }
}
