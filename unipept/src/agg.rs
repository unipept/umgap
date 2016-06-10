use taxon::{Taxon, TaxonId};

pub trait Aggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>) -> &Taxon;
}
