use taxon::{Taxon, TaxonId};

pub trait Aggregator {
    fn new(taxons: Vec<Taxon>, ranked_only: bool) -> Self;
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon;
}
