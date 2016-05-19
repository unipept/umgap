
use std::collections::HashMap;

use taxon;
use taxon::{TaxonId, Taxon};
use agg::Aggregator;

pub struct RTLCalculator {
    pub taxons: Vec<Option<Taxon>>,
    pub ancestors: Vec<Option<TaxonId>>,
}

impl Aggregator for RTLCalculator {
    fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self {
        // Views on the taxons
        let tree   = taxon::TaxonTree::new(&taxons);
        let mut by_id: Vec<Option<Taxon>> = (0..tree.max + 1).map(|_| None).collect();
        for taxon in taxons {
            let id = taxon.id;
            by_id[id] = Some(taxon);
        }

        // Precomputing
        let ancestors = if ranked_only {
            tree.filter_ancestors(|i: TaxonId| {
                let ref mtaxon = by_id[i];
                match *mtaxon {
                    None => false,
                    Some(ref taxon) => taxon.valid && taxon.rank != taxon::Rank::NoRank
                }
            })
        } else {
            tree.filter_ancestors(|i: TaxonId| {
                let ref mtaxon = by_id[i];
                match *mtaxon {
                    None => false,
                    Some(ref taxon) => taxon.valid
                }
            })
        };

        RTLCalculator {
            taxons: by_id,
            ancestors: ancestors,
        }

    }
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        // Count the occurences
        let mut counts = HashMap::new();
        for taxon in taxons { *counts.entry(*taxon).or_insert(0) += 1; }

        // current method: for each taxon id, loop to the root and add if the ancestor has a count
        // alternative method: for each taxon id, loop over all other and add if the other is an ancestor

        let mut rtl_counts = counts.clone();
        for (taxon, count) in rtl_counts.iter_mut() {
            let mut next = *taxon;
            while let Some(ancestor) = self.ancestors[next] {
                if ancestor == next { break; }
                *count += *counts.get(&next).unwrap_or(&0);
                next = ancestor;
            }
        }

        let rtl_taxon_id = *rtl_counts.iter()
                                      .max_by_key(|&(_, count)| count)
                                      .unwrap_or((&1, &0))
                                      .0;
        self.taxons[rtl_taxon_id].as_ref().expect("Taxonomy inconsistency.")
    }
}
