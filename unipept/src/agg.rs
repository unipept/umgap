//! Defines aggregation operations over a taxon tree.

use std::collections::HashMap;
use std::fmt;
use std::error;

use taxon::TaxonId;

/// Allows to aggregate over a taxon tree.
pub trait Aggregator {
    /// Aggregates a set of scored taxons into a resulting taxon id.
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> Result<TaxonId, Error>;

    /// Aggregates a list of taxons into a resulting taxon id.
    fn counting_aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, Error> {
        let taxons = taxons.iter().map(|&t| (t, 1.0));
        self.aggregate(&count(taxons))
    }
}

/// Returns how many times each taxon occurs in a vector of taxons.
pub fn count<T>(taxons: T) -> HashMap<TaxonId, f32>
        where T: Iterator<Item=(TaxonId, f32)> {
    let mut counts = HashMap::new();
    for (taxon, count) in taxons {
        *counts.entry(taxon).or_insert(0.0) += count;
    }
    counts
}

/// Defines an error which occurred during an aggregation operation
#[derive(Debug, PartialEq)]
pub enum Error {
    /// Returned if the collection to be aggregated is empty
    EmptyInput,
    /// Returned if the collection contains an unknown TaxonId
    UnknownTaxon(TaxonId),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::EmptyInput      => write!(f, "Aggregation called on empty list"),
            Error::UnknownTaxon(t) => write!(f, "Unknown Taxon ID: {}", t),
        }
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        match *self {
            Error::EmptyInput      => "The `aggregate` method was called on an empty list.",
            Error::UnknownTaxon(_) => "Result of aggregating taxa no in the list of known taxa.",
        }
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

impl From<TaxonId> for Error {
    fn from(id: TaxonId) -> Error {
        Error::UnknownTaxon(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use taxon::TaxonList;
    use fixtures;
    use rmq;
    use tree;

    fn aggregators(by_id: &TaxonList) -> Vec<Box<Aggregator>> { vec![
        Box::new(rmq::lca::LCACalculator::new(fixtures::tree())),
        Box::new(rmq::rtl::RTLCalculator::new(fixtures::ROOT, by_id)),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 0.0)),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 1.0)),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), 0.5)),
        Box::new(tree::lca::LCACalculator::new(fixtures::ROOT, by_id)),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 0.0)),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 1.0)),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, 0.5)),
    ] }

    #[test]
    fn test_empty_query() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_eq!(
                Err(Error::EmptyInput),
                aggregator.counting_aggregate(&Vec::new())
            );
        }
    }

    #[test]
    fn test_singleton_is_singleton() {
        for aggregator in aggregators(&fixtures::by_id()) {
            for taxon in fixtures::taxon_list() {
                assert_eq!(Ok(taxon.id), aggregator.counting_aggregate(&vec![taxon.id]));
            }
        }
    }

    #[test]
    fn test_invalid_taxa() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_eq!(
                Err(Error::UnknownTaxon(5)),
                aggregator.counting_aggregate(&vec![5])
            );
            assert_eq!(
                Err(Error::UnknownTaxon(5)),
                aggregator.counting_aggregate(&vec![1, 2, 5, 1])
            );
        }
    }
}
