//! Defines the possible errors in Unipept.

use std::io;
use std::num;

use csv;

use fst;

use regex;

use crate::agg;
use crate::dna::translation;
use crate::taxon;

/// Anything that can go wrong in umgap.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Error propagated from Taxon
    #[error(transparent)]
    Taxon(#[from] taxon::Error),
    /// Error propagated from Translation
    #[error(transparent)]
    Translation(#[from] translation::Error),
    /// Error propagated from Aggregation
    #[error(transparent)]
    Aggregation(#[from] agg::Error),
    /// CSV
    #[error(transparent)]
    Csv(#[from] csv::Error),
    /// Fst
    #[error(transparent)]
    Fst(#[from] fst::Error),
    /// IO
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Parse Integer
    #[error(transparent)]
    ParseI(#[from] num::ParseIntError),
    /// Parse Float
    #[error(transparent)]
    ParseF(#[from] num::ParseFloatError),
    /// Regex
    #[error(transparent)]
    Regex(#[from] regex::Error),
    /// Invalid invocation
    #[error("Invalid invocation: {0}")]
    InvalidInvocation(String),
    /// Anything else that went wrong
    #[error("{0}")]
    Msg(String),
}

// error-chain generated these, and a fair amount of code leans on them by doing
// `.map_err(|e| e.to_string())?` or `.ok_or("...")?`.
impl From<String> for Error {
    fn from(message: String) -> Self {
        Error::Msg(message)
    }
}

impl From<&str> for Error {
    fn from(message: &str) -> Self {
        Error::Msg(message.to_string())
    }
}

/// A Result with an umgap Error
pub type Result<T> = std::result::Result<T, Error>;
