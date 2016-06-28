
use std::collections::HashMap;
use std::fmt;
use std::error;

use taxon::TaxonId;

pub trait Aggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, Error>;
}

pub fn count(taxons: &Vec<TaxonId>) -> HashMap<TaxonId, usize> {
    let mut counts = HashMap::new();
    for taxon in taxons {
        *counts.entry(*taxon).or_insert(0) += 1;
    }
    counts
}

#[derive(Debug, PartialEq)]
pub enum Error {
    EmptyInput,
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
    extern crate num_rational;
    use self::num_rational::Ratio;
    use super::*;
    use taxon::Taxon;
    use fixtures;
    use rmq;
    use tree;

    fn aggregators(by_id: &Vec<Option<Taxon>>) -> Vec<Box<Aggregator>> { vec![
        Box::new(rmq::lca::LCACalculator::new(fixtures::tree())),
        Box::new(rmq::rtl::RTLCalculator::new(fixtures::ROOT, by_id)),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(0, 1))),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(1, 1))),
        Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(1, 2))),
        Box::new(tree::lca::LCACalculator::new(fixtures::ROOT, by_id)),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, Ratio::new(0, 1))),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, Ratio::new(1, 1))),
        Box::new(tree::mix::MixCalculator::new(fixtures::ROOT, by_id, Ratio::new(1, 2))),
    ] }

    #[test]
    fn test_empty_query() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_eq!(
                Err(Error::EmptyInput),
                aggregator.aggregate(&Vec::new())
            );
        }
    }

    #[test]
    fn test_singleton_is_singleton() {
        for aggregator in aggregators(&fixtures::by_id()) {
            for taxon in fixtures::taxon_list() {
                assert_eq!(Ok(taxon.id), aggregator.aggregate(&vec![taxon.id]));
            }
        }
    }

    #[test]
    fn test_invalid_taxa() {
        for aggregator in aggregators(&fixtures::by_id()) {
            assert_eq!(
                Err(Error::UnknownTaxon(5)),
                aggregator.aggregate(&vec![5])
            );
            assert_eq!(
                Err(Error::UnknownTaxon(5)),
                aggregator.aggregate(&vec![1, 2, 5, 1])
            );
        }
    }
}
