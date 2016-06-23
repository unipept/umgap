
pub mod taxon;
pub mod agg;
pub mod rmq;
pub mod tree;

#[cfg(test)]
pub mod fixtures;

pub const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
pub const PKG_NAME:    &'static str = env!("CARGO_PKG_NAME");
pub const PKG_AUTHORS: &'static str = env!("CARGO_PKG_AUTHORS");

