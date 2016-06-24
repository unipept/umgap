
use std::collections::HashMap;

use taxon::TaxonId;

pub trait Aggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, String>;
}

pub fn count(taxons: &Vec<TaxonId>) -> HashMap<TaxonId, usize> {
    let mut counts = HashMap::new();
    for taxon in taxons {
        *counts.entry(*taxon).or_insert(0) += 1;
    }
    counts
}

#[cfg(test)]
mod tests {
    extern crate num_rational;
    use self::num_rational::Ratio;
    use super::*;
    use fixtures;
    use rmq;
    use tree;


    #[test]
    fn test_empty_query() {
        let by_id = fixtures::by_id();
        let aggregators: Vec<Box<Aggregator>> = vec![
            Box::new(rmq::lca::LCACalculator::new(fixtures::tree())),
            Box::new(rmq::rtl::RTLCalculator::new(&by_id)),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(0, 1))),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(1, 1))),
            Box::new(rmq::mix::MixCalculator::new(fixtures::tree(), Ratio::new(1, 2))),
            Box::new(tree::lca::LCACalculator::new(fixtures::tree().root, &by_id)),
            Box::new(tree::mix::MixCalculator::new(fixtures::tree().root, &by_id, Ratio::new(0, 1))),
            Box::new(tree::mix::MixCalculator::new(fixtures::tree().root, &by_id, Ratio::new(1, 1))),
            Box::new(tree::mix::MixCalculator::new(fixtures::tree().root, &by_id, Ratio::new(1, 2))),
        ];
        for aggregator in aggregators {
            assert_eq!(
                Err("Aggregation called on empty list.".to_owned()),
                aggregator.aggregate(&Vec::new())
            );
        }
    }
}
