#![warn(missing_docs)]

pub mod taxon;
pub mod agg;
pub mod rmq;
pub mod tree;
pub mod errors;
pub mod io;

#[cfg(test)]
pub mod fixtures;

/// The current version
pub const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
/// The package name
pub const PKG_NAME:    &'static str = env!("CARGO_PKG_NAME");
/// The authors
pub const PKG_AUTHORS: &'static str = env!("CARGO_PKG_AUTHORS");

