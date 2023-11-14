//! Interface to pass a number of peptides (substrings of a protein/longer peptide)

/// A group of peptides
#[derive(Debug, PartialEq)]
pub struct Fragments {
    /// The header of a FASTA record
    pub header: String,
    /// The peptide from which fragments are taken
    pub peptide: String,
    /// The index and length of the fragments
    pub fragments: Vec<(usize, usize)>
}

impl Fragments {
    #[allow(missing_docs)]
    pub fn new(header: &str, peptide: &str, fragments: Vec<(usize, usize)>) -> Self {
        Fragments {
            header: header.to_string(),
            peptide: peptide.to_string(),
            fragments,
        }
    }
}
