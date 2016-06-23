extern crate num_rational;

use std::ops::Add;

use self::num_rational::Ratio;

use agg;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;
use tree::lca::LCACalculator;


pub struct MixCalculator {
    root: TaxonId,
    parents: Vec<Option<TaxonId>>,
    factor: Ratio<usize>
}

impl MixCalculator {
    pub fn new(root: TaxonId, taxonomy: &Vec<Option<Taxon>>, factor: Ratio<usize>) -> Self {
        let LCACalculator { root: r, parents: p } = LCACalculator::new(root, taxonomy);
        MixCalculator { factor: factor, root: r, parents: p }
    }
}

impl agg::Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> TaxonId {
        let counts  = agg::count(taxons);
        let subtree = SubTree::new(self.root, &self.parents, counts).unwrap()
            .collapse(&Add::add)
            .aggregate(&Add::add);

        let mut base = &subtree;
        while let Some(max) = base.children.iter().max_by_key(|c| c.value) {
            if Ratio::new(max.value, base.value) < self.factor { break; }
            base = max;
        }

        base.root
    }
}

