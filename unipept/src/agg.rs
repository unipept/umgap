
use std::collections::HashMap;

use taxon::{Taxon, TaxonId};

pub trait Aggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon;
}

pub fn count(taxons: &Vec<TaxonId>) -> HashMap<TaxonId, usize> {
    let mut counts = HashMap::new();
    for taxon in taxons {
        *counts.entry(*taxon).or_insert(0) += 1;
    }
    counts
}
