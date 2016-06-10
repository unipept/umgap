
// Serialize:
// https://doc.rust-lang.org/rustc-serialize/rustc_serialize/index.html

use std::collections::HashMap;

use rmq::RMQ;
use taxon::Taxon;
use taxon::TaxonId;
use taxon;
use agg::Aggregator;

pub struct LCACalculator {
    pub taxons: Vec<Option<Taxon>>,
    pub first_occurences: HashMap<TaxonId, usize>,
    pub euler_tour: Vec<TaxonId>,
    pub rmq_info: RMQ<usize>,
    pub snapping: Vec<Option<TaxonId>>,
}

impl LCACalculator {
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self {
        // Views on the taxons
        let length    = taxons.len();
        let tree      = taxon::TaxonTree::new(&taxons);
        let by_id     = taxon::taxa_vector_by_id(taxons);
        let snapping = tree.snapping(&by_id, ranked_only);

        // Euler tour
        let mut euler_tour       = Vec::with_capacity(length);
        let mut depths           = Vec::with_capacity(length);
        let mut first_occurences = HashMap::new();
        for (i, (taxon_id, depth)) in tree.into_iter().enumerate() {
            euler_tour.push(taxon_id);
            depths.push(depth);
            first_occurences.entry(taxon_id).or_insert(i);
        }

        // Result
        LCACalculator {
            taxons: by_id,
            first_occurences: first_occurences,
            euler_tour: euler_tour,
            rmq_info: RMQ::new(depths),
            snapping: snapping
        }
    }

    pub fn lca(&self, left: TaxonId, right: TaxonId) -> TaxonId {
        if left == right { return left; }
        let left_index  = *self.first_occurences.get(&left).expect("Unrecognized Taxon ID");
        let right_index = *self.first_occurences.get(&right).expect("Unrecognized Taxon ID");
        let rmq_index   =  self.rmq_info.query(left_index, right_index);
        self.snapping[self.euler_tour[rmq_index]].expect("LCA should be in the list.")
    }
}

impl Aggregator for LCACalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        let mut iter = taxons.into_iter().map(|&t| self.snapping[t].expect("Unrecognized Taxon ID"));
        let initial_level: Option<usize> = None;
        let first = iter.next().expect("LCA called on empty list");
        let (_, lca) = iter.fold((initial_level, first), |(join_level, left), right| {
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
        });
        let lca_taxon_id = self.snapping[lca].expect("Big bug: lca should exist.");
        self.taxons[lca_taxon_id].as_ref().expect("Taxonomy inconsistency.")
    }
}

