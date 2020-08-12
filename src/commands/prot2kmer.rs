//! The `umgap prot2kmer` command.

use std::borrow::Cow;
use std::io;

use crate::errors;
use crate::io::fasta;

/// The `umgap prot2kmer2lca` command takes one or more peptides as input and outputs all k-length
/// subsequences.
///
/// The input is given on *standard input* in a FASTA format. Per FASTA header should be a single
/// peptide, which may be hardwrapped with newlines. All overlapping k-mers of a peptide are written
/// to *standard output*, separated by newlines. (*k* is configurable via the `-k` option, and is 9
/// by default.)
///
///     $ cat input.fa
///     >header1
///     DAIGDVAKAYKKAG*S
///     $ umgap prot2kmer < input.fa
///     >header1
///     DAIGDVAKA
///     AIGDVAKAY
///     IGDVAKAYK
///     GDVAKAYKK
///     DVAKAYKKA
///     VAKAYKKAG
///     AKAYKKAG*
///     KAYKKAG*S
#[derive(Debug, StructOpt)]
pub struct ProtToKmer {
    /// The K in K-mers
    #[structopt(short = "k", long = "length", default_value = "9")]
    pub length: usize,
}

/// Implements the prot2kmer command.
pub fn prot2kmer(args: ProtToKmer) -> errors::Result<()> {
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), true).records() {
        let fasta::Record { header, sequence } = record?;
        if sequence[0].len() < args.length {
            continue;
        }
        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence[0]
                .as_bytes()
                .windows(args.length)
                .map(String::from_utf8_lossy)
                .map(Cow::into_owned)
                .collect(),
        })?;
    }
    Ok(())
}
