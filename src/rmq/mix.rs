//! Hybrid operation between MRL and LCA.

use std::collections::HashMap;
use std::collections::VecDeque;

extern crate ordered_float;
use self::ordered_float::NotNan;

use agg;
use rmq::lca::LCACalculator;
use taxon::{TaxonId, TaxonTree};

/// Struct capable of aggregating of a list of nodes in a TaxonTree, using a
/// hybrid approach between MRL and LCA. It can either prefer MRL or LCA more,
/// depending on a given ratio.
pub struct MixCalculator {
	lca_aggregator: LCACalculator,
	factor: f32,
}

#[derive(Clone, Copy)]
struct Weights {
	lca: f32,
	rtl: f32,
}

impl Weights {
	fn new() -> Self {
		Weights { lca: 0.0, rtl: 0.0 }
	}
}

impl MixCalculator {
	/// Constructs a new MixCalculator
	///
	/// # Arguments:
	/// * `taxonomy`, the TaxonTree containing all known taxons.
	/// * `factor`, a ratio (i.e. a number in [0.0, 1.0] which decides the ratio that MRL or LCA
	///   will be chosen as aggregation. If factor is 1, LCA will always be chosen; If factor is 0,
	///   MRL.
	pub fn new(taxonomy: TaxonTree, factor: f32) -> Self {
		MixCalculator { lca_aggregator: LCACalculator::new(taxonomy),
		                factor: factor }
	}
}

fn factorize(weights: Weights, factor: f32) -> f32 {
	weights.lca * factor + weights.rtl * (1.0 - factor)
}

impl agg::Aggregator for MixCalculator {
	fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> agg::Result<TaxonId> {
		let mut weights: HashMap<TaxonId, Weights> = HashMap::with_capacity(taxons.len());
		let mut queue: VecDeque<TaxonId> = taxons.keys().map(|&t| t).collect();

		// We collect all relevant taxa by starting out with all given taxa and iterate until no
		// more taxa are added. In each iteration, we add the lca of each pair already in the
		// collection.
		//
		// Each of these taxa is given two weights:
		// weight.lca is the count of all given taxa descending from this one (including itself).
		// weight.rtl is the count of all given taxa from which this one descends (including
		// itself).
		while let Some(left) = queue.pop_front() {
			if weights.contains_key(&left) {
				continue;
			}
			for (&right, &count) in taxons.iter() {
				let lca = self.lca_aggregator.lca(left, right)?;
				if lca == left || lca == right {
					let mut weight = weights.entry(left).or_insert(Weights::new());
					if lca == left {
						weight.lca += count;
					}
					if lca == right {
						weight.rtl += count;
					}
				} else {
					queue.push_back(lca);
				}
			}
		}

		weights.iter()
		       .max_by_key(|&(_, w)| NotNan::new(factorize(*w, self.factor)).unwrap())
		       .map(|tup| *tup.0)
		       .ok_or(agg::ErrorKind::EmptyInput.into())
	}
}

#[cfg(test)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod tests {
    use super::MixCalculator;
    use agg::Aggregator;
    use fixtures;

    #[test]
    fn test_full_rtl() {
        let aggregator = MixCalculator::new(fixtures::tree(), 0.0);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752, 185752]), Ok(185752));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]), Ok(10239));
    }

    #[test]
    fn test_full_lca() {
        let aggregator = MixCalculator::new(fixtures::tree(), 1.0);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]), Ok(1));
    }

    /* TODO: third example might fail because 12884 and 185751 have the same score. */
    #[test]
    fn test_one_half() {
        let aggregator = MixCalculator::new(fixtures::tree(), 0.5);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 12884, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 12884, 12884, 185751, 185752]), Ok(12884));
    }
}
