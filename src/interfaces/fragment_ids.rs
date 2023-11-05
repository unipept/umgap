//!
use std::fmt::Debug;

use crate::taxon::TaxonId;

/// A group of unaggregated fragment identifications
#[derive(Debug, PartialEq)]
pub struct FragmentIds {
    /// The header of the record
    pub header: String,
    /// The identifications of fragments
    pub taxa: Vec<TaxonId>,
}

impl FragmentIds {
    /// Constructor
    pub fn new(header: &str, taxa: Vec<TaxonId>) -> FragmentIds {
        FragmentIds {
            header: header.to_string(),
            taxa,
        }
    }
}
