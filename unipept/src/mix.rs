extern crate num_rational;
extern crate num_traits;

use std::collections::HashMap;
use std::collections::VecDeque;

use self::num_traits::identities::One;
use self::num_rational::Ratio;

use lca::LCACalculator;
use taxon::{Taxon, TaxonId};
use agg::Aggregator;

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
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool, factor: Ratio<usize>) -> Self {
        MixCalculator {
            lca_calculator: LCACalculator::new(taxons, ranked_only),
            factor: factor
        }
    }
}

fn factorize(weights: Weights, factor: Ratio<usize>) -> Ratio<usize> {
    Ratio::from_integer(weights.lca) * factor + Ratio::from_integer(weights.rtl) * (Ratio::one() - factor)
}

impl Aggregator for MixCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        // Count the occurences
        let mut counts = HashMap::new();
        for taxon in taxons {
            *counts.entry(*taxon).or_insert(0) += 1;
        }

        /* TODO: Collapsing of taxons with a single (present) child onto that child, to calculate
         * the LCA* in stead of the LCA. This calls for a tree implementation of the algorithm
         * below. */

        let mut weights: HashMap<TaxonId, Weights> = HashMap::with_capacity(counts.len());
        let mut queue:   VecDeque<TaxonId>         = counts.keys().map(|&t| t).collect();

        while let Some(left) = queue.pop_front() {
            if weights.contains_key(&left) { continue; }
            for (&right, &count) in counts.iter() {
                let lca = self.lca_calculator.lca(left, right);
                if lca == left || lca == right {
                    let mut weight = weights.entry(left).or_insert(Weights::new());
                    if lca == left  { weight.lca += count; println!("{} lies below {}, so adding {} to lca of {}", right, left, count, left); }
                    if lca == right { weight.rtl += count; println!("{} lies below {}, so adding {} to rtl of {}", left, right, count, left); }
                } else {
                    queue.push_back(lca);
                }
            }
        }

        for (taxon, weight) in weights.iter() {
            println!("{}: lca({}) and rtl({}), so {}", taxon, weight.lca, weight.rtl, factorize(*weight, self.factor));
        }

        let mix = weights.iter()
                         .max_by_key(|&(_, w)| factorize(*w, self.factor))
                         .expect("Should be at least one taxon")
                         .0;
        self.lca_calculator.taxons[*mix].as_ref().expect("None existing taxon.")
    }

}