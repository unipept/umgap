//! Defines aggregation operations over a taxon tree.

pub mod lineage;
pub mod rank;

use std::collections::HashMap;

use crate::taxon;
use crate::taxon::TaxonId;

/// Allows to aggregate over a taxon tree.
pub trait Aggregator {
    /// Aggregates a set of scored taxons into a resulting taxon id.
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> Result<TaxonId>;

    /// Aggregates a list of taxons into a resulting taxon id.
    fn counting_aggregate(&self, taxons: &[TaxonId]) -> Result<TaxonId> {
        let taxons = taxons.iter().map(|&t| (t, 1.0));
        self.aggregate(&count(taxons))
    }
}

/// Returns how many times each taxon occurs in a vector of taxons.
pub fn count<T>(taxons: T) -> HashMap<TaxonId, f32>
where
    T: Iterator<Item = (TaxonId, f32)>,
{
    let mut counts = HashMap::new();
    for (taxon, count) in taxons {
        *counts.entry(taxon).or_insert(0.0) += count;
    }
    counts
}

/// Filters any taxon in a frequency table with a frequency below the given amount.
pub fn filter(freq_table: HashMap<TaxonId, f32>, lower_bound: f32) -> HashMap<TaxonId, f32> {
    freq_table
        .into_iter()
        .filter(|&(_, freq)| freq >= lower_bound)
        .collect()
}

error_chain! {
    links {
        Taxon(taxon::Error, taxon::ErrorKind) #[doc = "Taxon"];
    }
    errors {
        /// Aggregation called on an empty list
        EmptyInput {
            description("Aggregration called on an empty list")
            display("Aggregration called on an empty list")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures;
    use crate::rmq;
    use crate::taxon::TaxonList;
    use crate::tree;

    fn aggregators(by_id: &TaxonList) -> Vec<Box<dyn Aggregator>> {
        vec![
            Box::new(rmq::lca::LCACalculator::new(fixtures::tree())),
            Box::new(rmq::rtl::RTLCalculator::new(fixtures::ROOT, by_id)),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 0.0)),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 1.0)),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 0.5)),
            Box::new(tree::lca::LCACalculator::new(fixtures::ROOT, by_id)),
            Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 0.0)),
            Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 1.0)),
            Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 0.5)),
        ]
    }

    #[test]
    fn test_empty_query() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_matches!(
                *aggregator
                    .counting_aggregate(&Vec::new())
                    .unwrap_err()
                    .kind(),
                ErrorKind::EmptyInput
            );
        }
    }

    #[test]
    fn test_singleton_is_singleton() {
        for aggregator in aggregators(&fixtures::by_id()) {
            for taxon in fixtures::taxon_list() {
                assert_matches!(aggregator.counting_aggregate(&vec![taxon.id]), Ok(tid) if tid == taxon.id);
            }
        }
    }

    #[test]
    fn test_invalid_taxa() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_matches!(
                *aggregator.counting_aggregate(&vec![5]).unwrap_err().kind(),
                ErrorKind::Taxon(taxon::ErrorKind::UnknownTaxon(5))
            );
            assert_matches!(
                *aggregator
                    .counting_aggregate(&vec![1, 2, 5, 1])
                    .unwrap_err()
                    .kind(),
                ErrorKind::Taxon(taxon::ErrorKind::UnknownTaxon(5))
            );
        }
    }
}
