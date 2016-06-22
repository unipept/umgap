extern crate num_rational;

use std::ops::Add;

use self::num_rational::Ratio;

use agg;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;
use tree::lca::LCACalculator;


pub struct MixCalculator {
    root: TaxonId,
    taxons: Vec<Option<Taxon>>,
    snapping: Vec<Option<TaxonId>>,
    factor: Ratio<usize>
}

impl MixCalculator {
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool, factor: Ratio<usize>) -> Self {
        let LCACalculator { root: r, taxons: t, snapping: s } = LCACalculator::new(taxons, ranked_only);
        MixCalculator { factor: factor, root: r, taxons: t, snapping: s }
    }
}

impl agg::Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        let counts  = agg::count(taxons);
        let subtree = SubTree::new(self.root, &self.taxons, counts).unwrap()
            .collapse(&Add::add)
            .aggregate(&Add::add);

        let mut base = &subtree;
        while let Some(max) = base.children.iter().max_by_key(|c| c.value) {
            if Ratio::new(max.value, base.value) < self.factor { break; }
            base = max;
        }

        let consensus = self.snapping[base.root].expect("LCA does not exist");
        self.taxons[consensus].as_ref().expect("No LCA* found.")
    }
}

