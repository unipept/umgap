//! The `umgap prot2tryp` command.

use std::io;

use crate::errors;
use crate::io::fasta;

#[structopt(verbatim_doc_comment)]
/// Splits the peptides in a FASTA stream into tryptic peptides
///
/// The `umgap prot2tryp` command takes one or more amino acid sequences as input and applies an *in
/// silico* trypsine digest.
///
/// The input is given in a FASTA format on *standard input* with a single peptide per FASTA header,
/// which may be hardwrapped with newlines. The peptides resulting from the digest are written
/// in FASTA format to *standard output*, with multiple peptides per FASTA header, separated by
/// newlines.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// AYKKAGVSGHVWQSDGITNCLLRGLTRVKEAVANRDSGNGYINKVYYWTVDKRATTRDALDAGVDGIMTNYPDVITDVLN
/// $ umgap prot2tryp tryptic-lca.index < input.fa
/// >header1
/// AYK
/// K
/// AGVSGHVWQSDGITNCLLR
/// GLTR
/// VK
/// EAVANR
/// DSGNGYINK
/// VYYWTVDK
/// R
/// ATTR
/// DALDAGVDGIMTNYPDVITDVLN
/// ```
///
/// Using the `-p` flag, you can change the splitting pattern. The default pattern `([KR])([^P])`
/// splits between any Lysine (K) or Arginine (R), followed by any amino acid that is not Proline
/// (P).
#[derive(Debug, StructOpt)]
pub struct ProtToTryp {
    /// The cleavage-pattern (regex), i.e. the pattern after which
    /// the next peptide will be cleaved for tryptic peptides)
    #[structopt(short = "p", long = "pattern", default_value = "([KR])([^P])")]
    pub pattern: String,
}

/// Implements the prot2tryp command
pub fn prot2tryp(args: ProtToTryp) -> errors::Result<()> {
    let pattern = regex::Regex::new(&args.pattern)?;

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), true).records() {
        let fasta::Record { header, sequence } = record?;

        // We will run the regex replacement twice, since a letter can be
        // matched twice (i.e. once after and once before the split).
        let first_run = pattern.replace_all(&sequence[0], "$1\n$2");

        writer.write_record(fasta::Record {
            header: header,
            sequence: pattern
                .replace_all(&first_run, "$1\n$2")
                .replace("*", "\n")
                .lines()
                .filter(|x| !x.is_empty())
                .map(ToOwned::to_owned)
                .collect(),
        })?;
    }
    Ok(())
}
