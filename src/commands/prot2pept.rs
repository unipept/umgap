//! The `umgap prot2pept` command.

use std::io;

use crate::errors;
use crate::io::fasta;

#[structopt(verbatim_doc_comment)]
/// The `umgap prot2pept` command takes one or more Amino Acid sequences as input, applies an *in
/// silico* peptide digest and outputs the result.
///
/// The input is given on *standard input*, in a FASTA format. Per FASTA header, it should contain
/// one sequence, possibly wrapped with newlines. It will write to *standard output*, in FASTA
/// format, the peptides resulting from the digest. Per FASTA header, there will be multiple
/// peptides, separated by newlines.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// AYKKAGVSGHVWQSDGITNCLLRGLTRVKEAVANRDSGNGYINKVYYWTVDKRATTRDALDAGVDGIMTNYPDVITDVLN
/// $ umgap prot2tryp2lca tryptic-lca.index < input.fa
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
/// Using the `-p` flag, you can change the splitting pattern (`([KR])([^P])`, so between a K and R
/// and before something that's not a P) to something else.
#[derive(Debug, StructOpt)]
pub struct ProtToPept {
    /// The cleavage-pattern (regex), i.e. the pattern after which
    /// the next peptide will be cleaved for tryptic peptides)
    #[structopt(short = "p", long = "pattern", default_value = "([KR])([^P])")]
    pub pattern: String,
}

/// Implements the prot2pept command
pub fn prot2pept(args: ProtToPept) -> errors::Result<()> {
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
