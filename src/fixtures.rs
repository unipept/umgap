use crate::rank::Rank;
use crate::taxon::*;

pub const ROOT: TaxonId = 1;
pub fn taxon_list() -> Vec<Taxon> {
    vec![
        Taxon::from_static(1, "root", Rank::NoRank, 1, true),
        Taxon::from_static(2, "Bacteria", Rank::Superkingdom, 1, true),
        Taxon::from_static(10239, "Viruses", Rank::Superkingdom, 1, true),
        Taxon::from_static(12884, "Viroids", Rank::Superkingdom, 1, true),
        Taxon::from_static(185751, "Pospiviroidae", Rank::Family, 12884, true),
        Taxon::from_static(185752, "Avsunviroidae", Rank::Family, 12884, true),
    ]
}

pub fn tree() -> TaxonTree {
    TaxonTree::new(&taxon_list())
}
pub fn by_id() -> TaxonList {
    TaxonList::new(taxon_list())
}
