
use std::ops::Add;

use agg::{Aggregator, count};
use taxon::{Taxon, TaxonId};
use tree::tree::SubTree;

pub struct LCACalculator {
    pub root: TaxonId,
    pub parents: Vec<Option<TaxonId>>,
}

impl LCACalculator {
    pub fn new(root: TaxonId, taxonomy: &Vec<Option<Taxon>>) -> Self {
        LCACalculator {
            root:    root,
            parents: taxonomy.iter().map(|mt| mt.as_ref().map(|t| t.parent)).collect()
        }
    }
}

impl Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> TaxonId {
        let counts  = count(taxons);
        let subtree = SubTree::new(self.root, &self.parents, counts).unwrap().collapse(&Add::add);
        subtree.root
    }
}

#[cfg(test)]
mod tests {
    use super::LCACalculator;
    use agg::Aggregator;
    use fixtures;

    #[test]
    fn test_two_on_same_path() {
        let calculator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(185752, calculator.aggregate(&vec![12884, 185752]));
        assert_eq!(185752, calculator.aggregate(&vec![185752, 12884]));
        assert_eq!(2, calculator.aggregate(&vec![1, 2]));
        assert_eq!(2, calculator.aggregate(&vec![2, 1]));
    }

    #[test]
    fn test_two_on_fork() {
        let calculator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(1, calculator.aggregate(&vec![2, 10239]));
        assert_eq!(1, calculator.aggregate(&vec![10239, 2]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 185751]));
    }

    #[test]
    fn test_three_on_triangle() {
        let calculator = LCACalculator::new(fixtures::tree().root, &fixtures::by_id());
        assert_eq!(12884, calculator.aggregate(&vec![12884, 185751, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![12884, 185752, 185751]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 12884, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 12884, 185751]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 185752, 12884]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 185751, 12884]));
    }
}
