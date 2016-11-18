//! Operations calculating the Maximal root-to-Leaf (MRL), i.e. the leaf that has most of
//! the given taxons on its path to root.

use taxon::{TaxonId, Taxon};
use agg;

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
    pub fn new(root: TaxonId, taxons: &Vec<Option<Taxon>>) -> Self {
        let mut ancestors: Vec<Option<TaxonId>> = taxons
            .iter()
            .map(|mtaxon| mtaxon.as_ref().map(|t| t.parent))
            .collect();
        ancestors[root] = None;

        RTLCalculator {
            root: root,
            ancestors: ancestors,
        }
    }
}

impl agg::Aggregator for RTLCalculator {
    /// Returns the taxon with the MRL for a given list of taxons.
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, agg::Error> {
        let counts         = agg::count(taxons);
        let mut rtl_counts = counts.clone();
        for (taxon, count) in rtl_counts.iter_mut() {
            let mut next = *taxon;
            while let Some(ancestor) = self.ancestors[next] {
                *count += *counts.get(&ancestor).unwrap_or(&0);
                next = ancestor;
            }
            if next != self.root {
                return Err(agg::Error::UnknownTaxon(next));
            }
        }

        rtl_counts.iter()
                  .max_by_key(|&(_, count)| count)
                  .map(|tup| *tup.0)
                  .ok_or(agg::Error::EmptyInput)
    }
}

#[cfg(test)]
mod tests {
    use super::RTLCalculator;
    use agg::Aggregator;
    use fixtures;


    #[test]
    fn test_all_on_same_path() {
        let calculator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_eq!(Ok(1), calculator.aggregate(&vec![1]));
        assert_eq!(Ok(12884), calculator.aggregate(&vec![1, 12884]));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![1, 12884, 185751]));
    }

    #[test]
    fn favouring_root() {
        let calculator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_eq!(Ok(185751), calculator.aggregate(&vec![1, 1, 1, 185751, 1, 1]));
    }

    #[test]
    fn leaning_close() {
        let calculator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert_eq!(Ok(185751), calculator.aggregate(&vec![1, 1, 185752, 185751, 185751, 1]));
    }

    #[test]
    fn non_deterministic() {
        let calculator = RTLCalculator::new(fixtures::ROOT, &fixtures::by_id());
        assert!(vec![Ok(185751), Ok(185752)].contains(&calculator.aggregate(&vec![1, 1, 185752, 185751, 1])))
    }
}
