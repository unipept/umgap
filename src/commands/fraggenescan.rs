//! The `umgap fraggenescan` command.

use std::path::PathBuf;

use crate::errors;

#[cfg_attr(rustfmt, rustfmt_skip)]
#[structopt(verbatim_doc_comment)]
/// Predicts genes found in a stream of DNA reads
///
/// The `umgap fraggenescan` command takes one or more DNA reads as input and predicts which fragments of
/// genes are encoded in the reads. By default it outputs the discovered protein fragments and their
/// location.
#[derive(Debug, StructOpt)]
pub struct FragGeneScan {
    /// File containing the training data of the prediction model
    #[structopt(short = "t", long = "training-file", parse(from_os_str))]
    pub training_file: PathBuf,

    /// The input reads are complete genomic sequences
    #[structopt(short = "w", long = "whole-genome")]
    pub whole_genome: bool,

    /// The number of threads to use
    #[structopt(short = "p", long = "threads", default_value = "1")]
    pub threads: usize,

    /// Write the untranslated DNA fragments to a file
    #[structopt(short = "d", long = "dna-file", parse(from_os_str))]
    pub dna_file: Option<PathBuf>,

    /// Write the prediction metadata of the prediction model to a file
    #[structopt(short = "e", long = "meta-file", parse(from_os_str))]
    pub meta_file: Option<PathBuf>,

}

/// Implements the fraggenescan command.
pub fn fraggenescan(args: FragGeneScan) -> errors::Result<()> {
    Ok(())
}
