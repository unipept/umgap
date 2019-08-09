//! Implementation of the UMGAP tools.

pub mod translate;

type Peptide = String;

/// A record, which is how data is passed between the UMGAP tools.
pub struct Record<T> {
	header: String,
	content: T
}
