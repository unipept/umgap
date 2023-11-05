//!

use crate::taxon::TaxonId;

/// A classification
pub struct Classification {
    /// The header of the record we're classifying
    pub header: String,
    /// The taxon the record is identified as
    pub taxon: TaxonId,
}
