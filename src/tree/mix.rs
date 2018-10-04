//! Hybrid operation between MRL and LCA.

use std::collections::HashMap;
use std::ops::Add;

extern crate ordered_float;
use self::ordered_float::NotNan;

use agg;
use taxon::{TaxonId, TaxonList};
use tree::lca::LCACalculator;
use tree::tree::Tree;


/// Struct capable of aggregating of a list of nodes in a TaxonTree, using a
/// hybrid approach between MRL and LCA. It can either prefer MRL or LCA more,
/// depending on a given ratio.
pub struct MixCalculator {
    root: TaxonId,
    parents: Vec<Option<TaxonId>>,
    factor: f32,
}

impl MixCalculator {
    /// Constructs a MixCalculator for a given taxon tree.
    ///
    /// # Arguments:
    /// * `root`     - the root of the taxon tree.
    /// * `taxonomy` - the taxons, indexed by their id.
    /// * `factor`   - A ratio (i.e. a number in [0.0, 1.0] which decides the
    ///                ratio that MRL or LCA will be chosen as aggregation.
    ///                If factor is 1, LCA will always be chosen; If factor is 0, MRL.
    pub fn new(root: TaxonId, taxonomy: &TaxonList, factor: f32) -> Self {
        let LCACalculator { root: r,
                            parents: p, } = LCACalculator::new(root, taxonomy);
        MixCalculator { factor: factor,
                        root: r,
                        parents: p, }
    }
}

impl agg::Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> agg::Result<TaxonId> {
        if taxons.len() == 0 {
            bail!(agg::ErrorKind::EmptyInput);
        }
        let subtree = Tree::new(self.root, &self.parents, taxons)?.collapse(&Add::add)
                                                                  .aggregate(&Add::add);

        let mut base = &subtree;
        while let Some(max) = base.children
                                  .iter()
                                  .max_by_key(|c| NotNan::new(c.value).unwrap())
        {
            if max.value / base.value < self.factor {
                break;
            }
            base = max;
        }

        Ok(base.root)
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
        let aggregator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), 0.0);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752, 185752]), Ok(185752));
        assert!(vec![185751, 185752].contains(&aggregator.counting_aggregate(&vec![1, 1, 10239, 10239, 12884, 185751, 185752]).unwrap()));
    }

    #[test]
    fn test_full_lca() {
        let aggregator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), 1.0);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]), Ok(1));
    }

    #[test]
    fn test_two_thirds() {
        let aggregator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), 0.66);
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 12884, 12884, 185751]), Ok(185751));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 12884, 10239, 185751, 185751, 185752]), Ok(12884));
    }
}
