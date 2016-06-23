
// Serialize:
// https://doc.rust-lang.org/rustc-serialize/rustc_serialize/index.html

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

