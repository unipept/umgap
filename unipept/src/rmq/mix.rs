extern crate num_rational;
extern crate num_traits;

use std::collections::HashMap;
use std::collections::VecDeque;

use self::num_traits::identities::One;
use self::num_rational::Ratio;

use rmq::lca::LCACalculator;
use taxon::{TaxonId, TaxonTree};
use agg::{Aggregator, count};

pub struct MixCalculator {
    lca_calculator: LCACalculator,
    factor: Ratio<usize>
}

#[derive(Clone, Copy)]
struct Weights {
    lca: usize,
    rtl: usize
}

impl Weights {
    fn new() -> Self {
        Weights { lca: 0, rtl: 0 }
    }
}

impl MixCalculator {
    pub fn new(taxonomy: TaxonTree, factor: Ratio<usize>) -> Self {
        MixCalculator {
            lca_calculator: LCACalculator::new(taxonomy),
            factor: factor
        }
    }
}

fn factorize(weights: Weights, factor: Ratio<usize>) -> Ratio<usize> {
    Ratio::from_integer(weights.lca) * factor + Ratio::from_integer(weights.rtl) * (Ratio::one() - factor)
}

impl Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> Result<TaxonId, String> {
        let     counts                             = count(taxons);
        let mut weights: HashMap<TaxonId, Weights> = HashMap::with_capacity(counts.len());
        let mut queue:   VecDeque<TaxonId>         = counts.keys().map(|&t| t).collect();

        while let Some(left) = queue.pop_front() {
            if weights.contains_key(&left) { continue; }
            for (&right, &count) in counts.iter() {
                let lca = self.lca_calculator.lca(left, right);
                if lca == left || lca == right {
                    let mut weight = weights.entry(left).or_insert(Weights::new());
                    if lca == left  { weight.lca += count; }
                    if lca == right { weight.rtl += count; }
                } else {
                    queue.push_back(lca);
                }
            }
        }

        weights.iter()
               .max_by_key(|&(_, w)| factorize(*w, self.factor))
               .map(|tup| *tup.0)
               .ok_or("Aggregation called on empty list.".to_owned())
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
        let calculator = MixCalculator::new(fixtures::tree(), Ratio::new(0, 1));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185752), calculator.aggregate(&vec![12884, 185751, 185752, 185752]));
        assert_eq!(Ok(10239), calculator.aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]));
    }

    #[test]
    fn test_full_lca() {
        let calculator = MixCalculator::new(fixtures::tree(), Ratio::new(1, 1));
        assert_eq!(Ok(12884), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(12884), calculator.aggregate(&vec![12884, 185751, 185752, 185752]));
        assert_eq!(Ok(1), calculator.aggregate(&vec![1, 1, 10239, 10239, 10239, 12884, 185751, 185752]));
    }

    /* TODO: third example might fail because 12884 and 185751 have the same score. */
    #[test]
    #[ignore]
    fn test_one_half() {
        let calculator = MixCalculator::new(fixtures::tree(), Ratio::new(1, 2));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![12884, 185751]));
        assert_eq!(Ok(185751), calculator.aggregate(&vec![1, 12884, 12284, 185751]));
    }
}
