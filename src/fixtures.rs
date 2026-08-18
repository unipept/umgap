//! A small taxonomy used to test the aggregation and tree operations.

use crate::rank::Rank;
use crate::taxon::*;

/// The taxon ID of the root of the test taxonomy.
pub const ROOT: TaxonId = 1;

/// The taxa making up the test taxonomy.
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

/// The test taxonomy as a tree.
pub fn tree() -> TaxonTree {
    TaxonTree::new(&taxon_list())
}

/// The test taxonomy indexed by taxon ID.
pub fn by_id() -> TaxonList {
    TaxonList::new(taxon_list())
}
