//! Allows calculating the Lowest Common Ancestor (LCA) using RMQ.

use std::collections::HashMap;

use agg;
use rmq::rmq::RMQ;
use taxon;
use taxon::TaxonId;

/// Struct capable of calculating the LCA of 2 nodes in a TaxonTree, using RMQ.
pub struct LCACalculator {
    /// Keeps track in which step of the euler tour a given taxon is encountered
    /// for the first time.
    pub first_occurences: HashMap<TaxonId, usize>,
    /// Contains the taxons in the order of an Eulerian tour.
    pub euler_tour: Vec<TaxonId>,
    /// The RMQ to which the LCA corresponds.
    pub rmq_info: RMQ<usize>,
}

impl LCACalculator {
    /// Creates an LCACalculator for the given TaxonTree.
    pub fn new(taxonomy: taxon::TaxonTree) -> Self {
        let mut euler_tour = Vec::new();
        let mut depths = Vec::new();
        let mut first_occurences = HashMap::new();
        for (i, (taxon_id, depth)) in taxonomy.into_iter().enumerate() {
            euler_tour.push(taxon_id);
            depths.push(depth);
            first_occurences.entry(taxon_id).or_insert(i);
        }

        // Result
        LCACalculator { first_occurences: first_occurences,
                        euler_tour: euler_tour,
                        rmq_info: RMQ::new(depths), }
    }

    /// Calculates the lowest common ancestor of 2 taxons.
    pub fn lca(&self, left: TaxonId, right: TaxonId) -> agg::Result<TaxonId> {
        let left_index = self.first_occurence(left)?;
        let right_index = self.first_occurence(right)?;
        let rmq_index = self.rmq_info.query(left_index, right_index);
        Ok(self.euler_tour[rmq_index])
    }

    fn first_occurence(&self, taxon_id: TaxonId) -> taxon::Result<usize> {
        self.first_occurences
            .get(&taxon_id)
            .ok_or(taxon::ErrorKind::UnknownTaxon(taxon_id).into())
            .map(|t| *t)
    }
}

impl agg::Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &HashMap<TaxonId, f32>) -> agg::Result<TaxonId> {
        if taxons.len() == 0 {
            bail!(agg::ErrorKind::EmptyInput);
        }
        let mut indices = taxons.keys().map(|t| self.first_occurence(*t));
        let mut consensus = indices.next().unwrap()?;
        let mut join_level = None::<usize>;
        for next_result in indices {
            let next = next_result?;
            if consensus == next {
                continue;
            }
            let rmq = self.rmq_info.query(consensus, next);
            let (mut lca, level) = match (rmq == consensus, rmq == next) {
                (false, false) => (rmq, Some(self.rmq_info.array[rmq])),
                (true, false) => (next, join_level),
                (false, true) => (consensus, join_level),
                (true, true) => panic!("Impossibru!"),
            };
            if join_level.map(|l| self.rmq_info.array[lca] > l)
                         .unwrap_or(false)
            {
                // join is below join level, we can't lower it
                lca = rmq;
            }
            consensus = lca;
            join_level = level
        }
        Ok(self.euler_tour[consensus])
    }
}

#[cfg(test)]
#[cfg_attr(rustfmt, rustfmt_skip)]
mod tests {
    use super::LCACalculator;
    use agg::Aggregator;
    use taxon::*;
    use fixtures;

    #[test]
    fn test_two_on_same_path() {
        let aggregator = LCACalculator::new(fixtures::tree());
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185752]), Ok(185752));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 12884]), Ok(185752));
        assert_matches!(aggregator.counting_aggregate(&vec![1, 2]), Ok(2));
        assert_matches!(aggregator.counting_aggregate(&vec![2, 1]), Ok(2));
    }

    #[test]
    fn test_two_on_fork() {
        let aggregator = LCACalculator::new(fixtures::tree());
        assert_matches!(aggregator.counting_aggregate(&vec![2, 10239]), Ok(1));
        assert_matches!(aggregator.counting_aggregate(&vec![10239, 2]), Ok(1));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 185751]), Ok(12884));
    }

    #[test]
    fn test_three_on_triangle() {
        let aggregator = LCACalculator::new(fixtures::tree());
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185751, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![12884, 185752, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 12884, 185752]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 12884, 185751]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185751, 185752, 12884]), Ok(12884));
        assert_matches!(aggregator.counting_aggregate(&vec![185752, 185751, 12884]), Ok(12884));
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
        let large_aggregator = LCACalculator::new(TaxonTree::new(&large_taxon_list()));
        assert_matches!(large_aggregator.counting_aggregate(&vec![9, 7]),  Ok(3));
        assert_matches!(large_aggregator.counting_aggregate(&vec![9, 10]), Ok(3));
        assert_matches!(large_aggregator.counting_aggregate(&vec![7, 9]),  Ok(3));
        assert_matches!(large_aggregator.counting_aggregate(&vec![14, 8]), Ok(3));
        assert_matches!(large_aggregator.counting_aggregate(&vec![14, 8]), Ok(3));
    }
}
