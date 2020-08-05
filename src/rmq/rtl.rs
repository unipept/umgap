//! Operations calculating the Maximal root-to-Leaf (MRL), i.e. the leaf that has most of
//! the given taxons on its path to root.

use std::collections::HashMap;

use self::ordered_float::NotNan;
use ordered_float;

use crate::agg;
use crate::taxon;
use crate::taxon::{TaxonId, TaxonList};

/// Struct capable of calculating the MRL of a list of nodes in a TaxonTree
pub struct RTLCalculator {
    /// The root node.
    pub root: TaxonId,
    /// Contains the parent for each node (to easily get from a node to the root).
    /// Nodes are indexed by their id.
    pub ancestors: Vec<Option<TaxonId>>,
}

impl RTLCalculator {
    /// Constructs an RTLCalculator.
    ///
    /// # Arguments:
    /// * `root`   - the root of the taxon tree
    /// * `taxons` - all taxons in the taxon tree, *indexed by their id*.
    pub fn new(root: TaxonId, taxons: &TaxonList) -> Self {
        let mut ancestors = taxons.ancestry();
        ancestors[root] = None;
        RTLCalculator {
            root: root,
            ancestors: ancestors,
        }
    }
}

impl agg::Aggregator for RTLCalculator {
    /// Returns the taxon with the MRL for a given list of taxons.
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> agg::Result<TaxonId> {
        let mut rtl_counts = taxons.clone();
        for (taxon, count) in rtl_counts.iter_mut() {
            let mut next = *taxon;
            while let Some(ancestor) = self.ancestors[next] {
                *count += *taxons.get(&ancestor).unwrap_or(&0.0);
                next = ancestor;
            }
            if next != self.root {
                bail!(agg::ErrorKind::Taxon(taxon::ErrorKind::UnknownTaxon(next)));
            }
        }

        rtl_counts
            .iter()
            .max_by_key(|&(_, &count)| NotNan::new(count).unwrap())
            .map(|tup| *tup.0)
            .ok_or(agg::ErrorKind::EmptyInput.into())
    }
}

#[cfg(test)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod tests {
    use super::RTLCalculator;
    use crate::agg::Aggregator;
    use crate::fixtures;


    #[test]
    fn test_all_on_same_path() {
        let aggregator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![1]), Ok(1));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 12884]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 12884, 185751]), Ok(185751));
    }

    #[test]
    fn favouring_root() {
        let aggregator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![1, 1, 1, 185751, 1, 1]), Ok(185751));
    }

    #[test]
    fn leaning_close() {
        let aggregator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![1, 1, 185752, 185751, 185751, 1]), Ok(185751));
    }

    #[test]
    fn non_deterministic() {
        let aggregator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert!(vec![185751, 185752].contains(&aggregator.counting_aggregate(&vec![1, 1, 185752, 185751, 1]).unwrap()));
    }
}
