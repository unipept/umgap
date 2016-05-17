
use taxon::TaxonId;

pub trait Aggregator {
    fn aggregate(&self, taxons: &Vec<TaxonId>, ranked_only: bool) -> TaxonId;
}

