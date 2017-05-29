//! Allows calculating the Lowest Common Ancestor (LCA).

use std::collections::HashMap;
use std::ops::Add;

use agg;
use taxon::{TaxonId, TaxonList};
use tree::tree::SubTree;

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
        LCACalculator {
            root:    root,
            parents: taxonomy.ancestry(),
        }
    }
}

impl agg::Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> Result<TaxonId, agg::Error> {
        if taxons.len() == 0 { bail!(agg::ErrorKind::EmptyInput); }
        let subtree = try!(SubTree::new(self.root, &self.parents, taxons)).collapse(&Add::add);
        Ok(subtree.root)
    }
}

#[cfg(test)]
mod tests {
    use super::LCACalculator;
    use agg::Aggregator;
    use fixtures;

    #[test]
    fn test_two_on_same_path() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(185752, aggregator.counting_aggregate(&vec![12884, 185752]).unwrap());
        assert_eq!(185752, aggregator.counting_aggregate(&vec![185752, 12884]).unwrap());
        assert_eq!(2, aggregator.counting_aggregate(&vec![1, 2]).unwrap());
        assert_eq!(2, aggregator.counting_aggregate(&vec![2, 1]).unwrap());
    }

    #[test]
    fn test_two_on_fork() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(1, aggregator.counting_aggregate(&vec![2, 10239]).unwrap());
        assert_eq!(1, aggregator.counting_aggregate(&vec![10239, 2]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185751, 185752]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185752, 185751]).unwrap());
    }

    #[test]
    fn test_three_on_triangle() {
        let aggregator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![12884, 185751, 185752]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![12884, 185752, 185751]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185751, 12884, 185752]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185752, 12884, 185751]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185751, 185752, 12884]).unwrap());
        assert_eq!(12884, aggregator.counting_aggregate(&vec![185752, 185751, 12884]).unwrap());
    }
}
