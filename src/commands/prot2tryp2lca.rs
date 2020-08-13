//! The `umgap prot2tryp2lca` command.

use std::collections::HashSet;
use std::fs;
use std::io;
use std::path::PathBuf;

use fst;

use regex;

use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::errors;
use crate::io::fasta;

#[structopt(verbatim_doc_comment)]
/// The `umgap prot2tryp2lca` command takes one or more peptides as input, splits those into
/// tryptic peptides, possibly filtering them, and outputs their lowest common ancestors. It is a
/// combination of the `umgap prot2pept`, `umgap filter` and `umgap pept2lca` commands to allow more
/// efficient parallel computing (c.f. their documentation for details).
///
/// The input is given on *standard input* in a FASTA format. Per FASTA header should be a single
/// peptide, which may be hardwrapped with newlines. The command prints the lowest common ancestors
/// for each tryptic peptide found in this peptide to standard output.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// AYKKAGVSGHVWQSDGITNCLLRGLTRVKEAVANRDSGNGYINKVYYWTVDKRATTRDALDAGVDGIMTNYPDVITDVLN
/// $ umgap prot2tryp2lca tryptic-lca.index < input.fa
/// >header1
/// 571525
/// 1
/// 571525
/// 6920
/// ```
#[derive(Debug, StructOpt)]
pub struct ProtToTrypToLca {
    /// Map unknown sequences to 0 instead of ignoring them
    #[structopt(short = "o", long = "one-on-one")]
    pub one_on_one: bool,

    /// An FST to query
    #[structopt(parse(from_os_str))]
    pub fst_file: PathBuf,

    /// Load FST in memory instead of mmap'ing the file contents. This makes
    /// querying significantly faster, but requires some time to load the FST
    /// into memory.
    #[structopt(short = "m", long = "in-memory")]
    pub fst_in_memory: bool,

    /// How much reads to group together in one chunk, bigger chunks decrease
    /// the overhead caused by multithreading. Because the output order is not
    /// necessarily the same as the input order, having a chunk size which is
    /// a multiple of 12 (all 6 translations multiplied by the two paired-end
    /// reads) will keep records of the same reads together.
    #[structopt(short = "c", long = "chunksize", default_value = "240")]
    pub chunk_size: usize,

    /// The cleavage-pattern (regex), i.e. the pattern after which
    /// the next peptide will be cleaved for tryptic peptides)
    #[structopt(short = "p", long = "pattern", default_value = "([KR])([^P])")]
    pub pattern: String,

    /// Minimum length of tryptic peptides to be mapped
    #[structopt(short = "l", long = "minlen", default_value = "5")]
    pub min_length: usize,

    /// Maximum length of tryptic peptides to be mapped
    #[structopt(short = "L", long = "maxlen", default_value = "50")]
    pub max_length: usize,

    /// Keep only tryptic peptides containing these letters
    #[structopt(short = "k", long = "keep", default_value = "")]
    pub contains: String,

    /// Drop tryptic peptides containing these letters
    #[structopt(short = "d", long = "drop", default_value = "")]
    pub lacks: String,
}

/// Implements the prot2tryp2lca command.
pub fn prot2tryp2lca(args: ProtToTrypToLca) -> errors::Result<()> {
    let fst = if args.fst_in_memory {
        let bytes = fs::read(&args.fst_file)?;
        fst::Map::from_bytes(bytes)?
    } else {
        unsafe { fst::Map::from_path(&args.fst_file) }?
    };
    let default = if args.one_on_one { Some(0) } else { None };
    let pattern = regex::Regex::new(&args.pattern)?;
    let contains = args.contains.chars().collect::<HashSet<char>>();
    let lacks = args.lacks.chars().collect::<HashSet<char>>();

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
                    // We will run the regex replacement twice, since a letter can be
                    // matched twice (i.e. once after and once before the split).
                    let first_run = pattern.replace_all(&seq, "$1\n$2");
                    for peptide in pattern
                        .replace_all(&first_run, "$1\n$2")
                        .replace("*", "\n")
                        .lines()
                        .filter(|x| !x.is_empty())
                        .filter(|seq| {
                            let length = seq.len();
                            length >= args.min_length && length <= args.max_length
                        })
                        .filter(|seq| {
                            (contains.is_empty() && lacks.is_empty()) || {
                                let set = seq.chars().collect::<HashSet<char>>();
                                contains.intersection(&set).count() == contains.len()
                                    && lacks.intersection(&set).count() == 0
                            }
                        })
                    {
                        if let Some(lca) = fst.get(&peptide).map(Some).unwrap_or(default) {
                            chunk_output.push_str(&format!("{}\n", lca));
                        }
                    }
                }
            }
            // TODO: make this the result of the map
            // and print using a Writer
            print!("{}", chunk_output);
            Ok(())
        })
        .collect()
}
