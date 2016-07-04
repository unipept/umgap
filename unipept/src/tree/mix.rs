
use std::ops::Add;

extern crate num_rational;
use self::num_rational::Ratio;

use agg;
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;
use tree::lca::LCACalculator;


pub struct MixCalculator {
    root: TaxonId,
    parents: Vec<Option<TaxonId>>,
    factor: Ratio<usize>
}

impl MixCalculator {
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

