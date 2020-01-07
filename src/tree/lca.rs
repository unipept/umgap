//! Allows calculating the Lowest Common Ancestor (LCA).

use std::collections::HashMap;
use std::ops::Add;

use crate::agg;
use crate::taxon::{TaxonId, TaxonList};
use crate::tree::tree::Tree;

/// Struct capable of calculating the LCA of 2 nodes in a TaxonTree.
pub struct LCACalculator {
	/// The root of the taxon tree.
	pub root: TaxonId,
	/// Contains the ancestor for each node. Nodes are indexed by their id.
	pub parents: Vec<Option<TaxonId>>,
}

impl LCACalculator {
	/// Constructs an LCACalculator for a given taxon tree.
	///
	/// # Arguments:
	/// * `root`     - the root of the taxon tree.
	/// * `taxonomy` - the taxons, indexed by their id.
	pub fn new(root: TaxonId, taxonomy: &TaxonList) -> Self {
		LCACalculator { root: root,
		                parents: taxonomy.ancestry() }
	}
}

impl agg::Aggregator for LCACalculator {
	fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> Result<TaxonId, agg::Error> {
		if taxons.len() == 0 {
			bail!(agg::ErrorKind::EmptyInput);
		}
		let subtree = Tree::new(self.root, &self.parents, taxons)?.collapse(&Add::add);
		Ok(subtree.root)
	}
}

#[cfg(test)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod tests {
    use super::LCACalculator;
    use crate::agg::Aggregator;
    use crate::fixtures;

    #[test]
    fn test_two_on_same_path() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185752]), Ok(185752));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 12884]), Ok(185752));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 2]), Ok(2));
        assert_matches!(aggregator.counting_aggregate(&vec![2, 1]), Ok(2));
    }

    #[test]
    fn test_two_on_fork() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![2, 10239]), Ok(1));
        assert_matches!(aggregator.counting_aggregate(&vec![10239, 2]), Ok(1));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 185751]), Ok(12884));
    }

    #[test]
    fn test_three_on_triangle() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185752, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 12884, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 12884, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 185752, 12884]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 185751, 12884]), Ok(12884));
    }
}
