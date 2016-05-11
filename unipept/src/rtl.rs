
use std::collections::HashMap;

use taxon::{TaxonId, Taxon};
use lca::LCACalculator;
use agg::Aggregator;

pub struct RTLCalculator {
    pub lca_calculator: LCACalculator,
}

impl RTLCalculator {
    pub fn new(taxons: Vec<Taxon>) -> RTLCalculator {
        RTLCalculator { lca_calculator: LCACalculator::new(taxons) }
    }
}

impl Aggregator for RTLCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>, ranked_only: bool) -> TaxonId {
        // Count the occurences
        let mut counts = HashMap::new();
        for taxon in taxons { *counts.entry(*taxon).or_insert(0) += 1; }

        // current method: for each taxon id, loop to the root and add if the ancestor has a count
        // alternative method: for each taxon id, loop over all other and add if the other is an ancestor
        let ancestors = if ranked_only {
            &self.lca_calculator.ranked_ancestors
        } else {
            &self.lca_calculator.valid_ancestors
        };

        let mut rtl_counts = counts.clone();
        for (taxon, count) in rtl_counts.iter_mut() {
            let mut next = *taxon;
            while let Some(ancestor) = ancestors[next] {
                if ancestor == next { break; }
                *count += *counts.get(&next).unwrap_or(&0);
                next = ancestor;
            }
        }

        *rtl_counts.iter()
                   .max_by_key(|&(_, count)| count)
                   .unwrap_or((&1, &0))
                   .0
    }
}
