//! Contains the implementation of several programs used by the
//! [Unipept pipeline](http://unipept.ugent.be/).

#![doc(
    html_logo_url = "http://unipept.ugent.be/logo.png",
    html_favicon_url = "http://unipept.ugent.be/favicon.ico"
)]
#![deny(missing_docs)]
#![recursion_limit = "128"]

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate lazy_static;

#[cfg(test)]
#[macro_use]
extern crate assert_matches;

#[macro_use]
extern crate structopt;

#[macro_use]
extern crate strum_macros;

pub mod agg;
pub mod cli;
pub mod dna;
pub mod errors;
pub mod io;
pub mod rank;
pub mod rmq;
pub mod taxon;
pub mod tree;
pub mod utils;

#[cfg(test)]
pub mod fixtures;

/// The current version
pub const PKG_VERSION: &str = env!("CARGO_PKG_VERSION");
/// The package name
pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");
/// The authors
pub const PKG_AUTHORS: &str = env!("CARGO_PKG_AUTHORS");
