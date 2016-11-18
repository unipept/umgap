//! Hybrid operation between MRL and LCA.

use std::ops::Add;

extern crate num_rational;
use self::num_rational::Ratio;

use agg;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;
use tree::lca::LCACalculator;


/// Struct capable of aggregating of a list of nodes in a TaxonTree, using a
/// hybrid approach between MRL and LCA. It can either prefer MRL or LCA more,
/// depending on a given ratio.
pub struct MixCalculator {
    root: TaxonId,
    parents: Vec<Option<TaxonId>>,
    factor: Ratio<usize>
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
    pub fn new(root: TaxonId, taxonomy: &Vec<Option<Taxon>>, factor: Ratio<usize>) -> Self {
        let LCACalculator { root: r, parents: p } = LCACalculator::new(root, taxonomy);
        MixCalculator { factor: factor, root: r, parents: p }
    }
}

impl agg::Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, agg::Error> {
        if taxons.len() == 0 { return Err(agg::Error::EmptyInput); }
        let counts  = agg::count(taxons);
        let subtree = try!(SubTree::new(self.root, &self.parents, counts))
            .collapse(&Add::add)
            .aggregate(&Add::add);

        let mut base = &subtree;
        while let Some(max) = base.children.iter().max_by_key(|c| c.value) {
            println!("{} <? {}", Ratio::new(max.value, base.value), self.factor);
            if Ratio::new(max.value, base.value) < self.factor { break; }
            base = max;
        }

        Ok(base.root)
    }
}

#[cfg(test)]
mod tests {
    extern crate num_rational;
    use self::num_rational::Ratio;

    use super::MixCalculator;
    use agg::Aggregator;
    use fixtures;

    #[test]
    fn test_full_rtl() {
        let calculator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), Ratio::new(0, 1));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185752), calculator.aggregate(&vec![12884, 185751, 185752, 185752]));
        assert!(vec![Ok(185751), Ok(185752)].contains(&calculator.aggregate(&vec![1, 1, 10239, 10239, 12884, 185751, 185752])));
    }

    #[test]
    fn test_full_lca() {
        let calculator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), Ratio::new(1, 1));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(12884), calculator.aggregate(&vec![12884, 185751, 185752, 185752]));
        assert_eq!(Ok(1), calculator.aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]));
    }

    #[test]
    fn test_two_thirds() {
        let calculator = MixCalculator::new(fixtures::ROOT, &fixtures::by_id(), Ratio::new(2, 3));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![1, 12884, 12884, 185751]));
        assert_eq!(Ok(12884), calculator.aggregate(&vec![1, 12884, 10239, 185751, 185751, 185752]));
    }
}

