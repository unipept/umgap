//! The `umgap pept2lca` command.

use std::fs;
use std::io;
use std::path::PathBuf;

use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::errors;
use crate::io::fasta;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Maps a FASTA stream of peptides to taxon IDs
///
/// The `umgap pept2lca` command takes one or more amino acid sequences and looks up the
/// corresponding taxon ID in an index file (as build by the `umgap buildindex` command).
///
/// The input is given in FASTA format on *standard input*. Per FASTA header, there can be multiple
/// sequences, each on a line. In the following example we match tryptic peptides on their lowest
/// common ancestor in the NCBI taxonomy.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// AAALTER
/// ENFVYLAK
/// $ umgap pept2lca tryptic-peptides.index < input.fa
/// >header1
/// 2
/// 3398
/// ```
///
/// By default, sequences not found in the index are ignored. Using the `-o` (`--on-on-one`) flag,
/// they are mapped to 0, instead.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// NOTATRYPTICPEPTIDE
/// ENFVYLAK
/// $ umgap pept2lca -o tryptic-peptides.index < input.fa
/// >header1
/// 0
/// 3398
/// ```
pub struct PeptToLca {
    /// Map unknown sequences to 0 instead of ignoring them
    #[structopt(short = "o", long = "one-on-one")]
    pub one_on_one: bool,

    /// An index that maps peptides to taxon IDs
    #[structopt(parse(from_os_str))]
    pub fst_file: PathBuf,

    /// Load index in memory instead of memory mapping the file contents. This
    /// makes querying significantly faster, but requires some initialization
    /// time.
    #[structopt(short = "m", long = "in-memory")]
    pub fst_in_memory: bool,

    /// Number of reads grouped into one chunk. Bigger chunks decrease
    /// the overhead caused by multithreading. Because the output order is not
    /// necessarily the same as the input order, having a chunk size which is
    /// a multiple of 12 (all 6 translations multiplied by the two paired-end
    /// reads) will keep FASTA records that originate from the same reads
    /// together.
    #[structopt(short = "c", long = "chunksize", default_value = "240")]
    pub chunk_size: usize,
}

/// Implements the pept2lca command
pub fn pept2lca(args: PeptToLca) -> errors::Result<()> {
    let fst = if args.fst_in_memory {
        let bytes = fs::read(args.fst_file)?;
        fst::Map::from_bytes(bytes)?
    } else {
        unsafe { fst::Map::from_path(args.fst_file) }?
    };

    let default = if args.one_on_one { Some(0) } else { None };

    fasta::Reader::new(io::stdin(), false)
        .records()
        .chunked(args.chunk_size)
        .par_bridge()
        .map(|chunk| {
            let chunk = chunk?;
            let mut chunk_output = String::new();
            for read in chunk {
                chunk_output.push_str(&format!(">{}\n", read.header));
                for seq in read.sequence {
                    if let Some(lca) = fst.get(&seq).map(Some).unwrap_or(default) {
                        chunk_output.push_str(&format!("{}\n", lca));
                    }
                }
            }
            print!("{}", chunk_output);
            Ok(())
        })
        .collect()
}
