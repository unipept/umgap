//! Contains the implementation of several programs used by the
//! [Unipept pipeline](http://unipept.ugent.be/).

#![doc(html_logo_url = "http://unipept.ugent.be/logo.png",
       html_favicon_url = "http://unipept.ugent.be/favicon.ico")]
#![deny(missing_docs)]

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate lazy_static;

#[cfg(test)]
#[macro_use]
extern crate assert_matches;

extern crate regex;

pub mod taxon;
pub mod agg;
pub mod rmq;
pub mod tree;
pub mod errors;
pub mod io;
pub mod dna;

#[cfg(test)]
pub mod fixtures;

/// The current version
pub const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
/// The package name
pub const PKG_NAME:    &'static str = env!("CARGO_PKG_NAME");
/// The authors
pub const PKG_AUTHORS: &'static str = env!("CARGO_PKG_AUTHORS");

