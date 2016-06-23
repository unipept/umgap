
use std::collections::HashMap;

use rmq::rmq::RMQ;
use taxon;
use taxon::TaxonId;
use agg::Aggregator;

pub struct LCACalculator {
    pub first_occurences: HashMap<TaxonId, usize>,
    pub euler_tour: Vec<TaxonId>,
    pub rmq_info: RMQ<usize>,
}

impl LCACalculator {
    pub fn new(taxonomy: taxon::TaxonTree) -> Self {
        let mut euler_tour       = Vec::new();
        let mut depths           = Vec::new();
        let mut first_occurences = HashMap::new();
        for (i, (taxon_id, depth)) in taxonomy.into_iter().enumerate() {
            euler_tour.push(taxon_id);
            depths.push(depth);
            first_occurences.entry(taxon_id).or_insert(i);
        }

        // Result
        LCACalculator {
            first_occurences: first_occurences,
            euler_tour: euler_tour,
            rmq_info: RMQ::new(depths),
        }
    }

    pub fn lca(&self, left: TaxonId, right: TaxonId) -> TaxonId {
        if left == right { return left; }
        let left_index  = *self.first_occurences.get(&left).expect("Unrecognized Taxon ID");
        let right_index = *self.first_occurences.get(&right).expect("Unrecognized Taxon ID");
        let rmq_index   =  self.rmq_info.query(left_index, right_index);
        self.euler_tour[rmq_index]
    }
}

impl Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> TaxonId {
        let mut iter = taxons.into_iter();
        let initial_level: Option<usize> = None;
        let first = *iter.next().expect("LCA called on empty list");
        iter.fold((initial_level, first), |(join_level, left), &right| {
            if left == right { return (join_level, left); }
            let left_index  = *self.first_occurences.get(&left).expect("Unrecognized Taxon ID");
            let right_index = *self.first_occurences.get(&right).expect("Unrecognized Taxon ID");
            let rmq_index   = self.rmq_info.query(left_index, right_index);
            let (mut lca_index, level) = match (rmq_index == left_index, rmq_index == right_index) {
                (false, false) => (rmq_index, Some(self.rmq_info.array[rmq_index])),
                (true,  false) => (right_index, join_level),
                (false, true)  => (left_index, join_level),
                (true,  true)  => panic!("Impossibru!")
            };
            if join_level.map(|l| self.rmq_info.array[lca_index] > l).unwrap_or(false) {
                lca_index = rmq_index;
            }
            (level, self.euler_tour[lca_index])
        }).1
    }
}

#[cfg(test)]
mod tests {
    use super::LCACalculator;
    use agg::Aggregator;
    use taxon::*;
    use fixtures;

    #[test]
    fn test_two_on_same_path() {
        let calculator = LCACalculator::new(fixtures::tree());
        assert_eq!(185752, calculator.aggregate(&vec![12884, 185752]));
        assert_eq!(185752, calculator.aggregate(&vec![185752, 12884]));
        assert_eq!(2, calculator.aggregate(&vec![1, 2]));
        assert_eq!(2, calculator.aggregate(&vec![2, 1]));
    }

    #[test]
    fn test_two_on_fork() {
        let calculator = LCACalculator::new(fixtures::tree());
        assert_eq!(1, calculator.aggregate(&vec![2, 10239]));
        assert_eq!(1, calculator.aggregate(&vec![10239, 2]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 185751]));
    }

    #[test]
    fn test_three_on_triangle() {
        let calculator = LCACalculator::new(fixtures::tree());
        assert_eq!(12884, calculator.aggregate(&vec![12884, 185751, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![12884, 185752, 185751]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 12884, 185752]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 12884, 185751]));
        assert_eq!(12884, calculator.aggregate(&vec![185751, 185752, 12884]));
        assert_eq!(12884, calculator.aggregate(&vec![185752, 185751, 12884]));
    }

    fn taxon(id: TaxonId, parent: TaxonId) -> Taxon {
        Taxon::from_static(id, "", Rank::NoRank, parent, true)
    }

    fn large_taxon_list() -> Vec<Taxon> {
        vec![
            taxon(1, 1),
              taxon(2, 1),
                taxon(5, 2),
                taxon(6, 2),
              taxon(3, 1),
                taxon(7, 3),
                  taxon(10, 7),
                    taxon(13, 10),
                      taxon(14, 13),
                  taxon(15, 3),
                taxon(8, 3),
                  taxon(11, 8),
                  taxon(12, 8),
                taxon(9, 3),
              taxon(4, 1)
        ]
    }

    #[test]
    fn test_with_deeper_interns() {
        let large_calculator = LCACalculator::new(TaxonTree::new(&large_taxon_list()));
        assert_eq!(3, large_calculator.aggregate(&vec![9, 7]));
        assert_eq!(3, large_calculator.aggregate(&vec![9, 10]));
        assert_eq!(3, large_calculator.aggregate(&vec![7, 9]));
        assert_eq!(3, large_calculator.aggregate(&vec![14, 8]));
        assert_eq!(3, large_calculator.aggregate(&vec![14, 8]));
    }
}

