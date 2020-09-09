//! The `umgap filter` command.

use std::collections::HashSet;
use std::io;

use crate::errors;
use crate::io::fasta;

#[structopt(verbatim_doc_comment)]
/// Filters a FASTA stream of peptides
///
/// The `umgap filter` command takes a FASTA stream of peptides as input and outputs a filtered
/// stream.
///
/// The input is given in FASTA format on *standard input*. Per FASTA header, there may be multiple
/// peptides separated by newlines. Each of these peptides are checked against the filter criteria
/// and written to *standard output* they pass them. The criteria are specified as options:
///
/// * `-m 5` sets the minimum length of the peptides to 5 (which is the default).
/// * `-M 50` sets the maximum length of the peptides to 50 (which is the default).
/// * `-c LIK` requires the peptides to contain amino acids Leucine (L), Isoleucine (I) and Lysine
///   (K) (none by default).
/// * `-l LIK` removes the peptides containing the amino acids Leucine, Isoleucine or Lysine (none
///   by default).
///
/// ```sh
/// $ cat input.fa
/// >header1
/// AYKKAGVSGHVWQSDGITNCLLRGLTRVKEAVANRDSGNGYINKVYYWTVDKRATTRDALDAGVDGIMTNYPDVITDVLN
/// AYK
/// K
/// AGVSGHVWQSDGITNCLLR
/// GLTR
/// VK
/// EAVANR
/// DSGNGYINK
/// $ umgap filter < input.fa
/// >header1
/// AGVSGHVWQSDGITNCLLR
/// EAVANR
/// DSGNGYINK
/// $ umgap filter -m 0 -c R -l K < input.fa
/// >header1
/// AGVSGHVWQSDGITNCLLR
/// GLTR
/// EAVANR
/// ```
#[derive(Debug, StructOpt)]
pub struct Filter {
    /// Minimum length required
    #[structopt(short = "m", long = "minlen", default_value = "5")]
    pub min_length: usize,

    /// Maximum length allowed
    #[structopt(short = "M", long = "maxlen", default_value = "50")]
    pub max_length: usize,

    /// Amino acid symbols that a sequence must contain
    #[structopt(short = "c", long = "contains", default_value = "")]
    pub contains: String,

    /// Amino acid symbols that a sequence may not contain
    #[structopt(short = "l", long = "lacks", default_value = "")]
    pub lacks: String,
}

/// Implements the filter command.
pub fn filter(args: Filter) -> errors::Result<()> {
    let contains = args.contains.chars().collect::<HashSet<char>>();
    let lacks = args.lacks.chars().collect::<HashSet<char>>();

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let fasta::Record { header, sequence } = record?;

        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence
                .into_iter()
                .filter(|seq| {
                    let length = seq.len();
                    length >= args.min_length && length <= args.max_length
                })
                .filter(|seq| {
                    let set = seq.chars().collect::<HashSet<char>>();
                    contains.intersection(&set).count() == contains.len()
                        && lacks.intersection(&set).count() == 0
                })
                .collect(),
        })?;
    }
    Ok(())
}
