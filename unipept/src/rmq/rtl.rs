
use taxon;
use taxon::{TaxonId, Taxon};
use agg::{Aggregator, count};

pub struct RTLCalculator {
    pub taxons: Vec<Option<Taxon>>,
    pub ancestors: Vec<Option<TaxonId>>,
}

impl RTLCalculator {
    pub fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self {
        // Views on the taxons
        let tree      = taxon::TaxonTree::new(&taxons);
        let by_id     = taxon::taxa_vector_by_id(taxons);
        let snapping  = tree.snapping(&by_id, ranked_only);

        let mut ancestors: Vec<Option<TaxonId>> = by_id
            .iter()
            .map(|mtaxon| mtaxon.as_ref().and_then(|t| snapping[t.parent]))
            .collect();
        ancestors[1] = None;

        RTLCalculator {
            taxons: by_id,
            ancestors: ancestors
        }
    }
}

impl Aggregator for RTLCalculator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon {
        let counts         = count(taxons);
        let mut rtl_counts = counts.clone();
        for (taxon, count) in rtl_counts.iter_mut() {
            let mut next = *taxon;
            while let Some(ancestor) = self.ancestors[next] {
                *count += *counts.get(&ancestor).unwrap_or(&0);
                next = ancestor;
            }
        }

        let rtl_taxon_id = *rtl_counts.iter()
                                      .max_by_key(|&(_, count)| count)
                                      .unwrap_or((&1, &0))
                                      .0;
        self.taxons[rtl_taxon_id].as_ref().expect("Taxonomic inconsistency.")
    }
}
