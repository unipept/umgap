//! The `umgap prot2kmer` command.

use std::borrow::Cow;
use std::io;

use crate::errors;
use crate::io::fasta;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Splits a FASTA stream of peptides into k-mers
///
/// The `umgap prot2kmer` command takes one or more peptides as input and outputs all their k-length
/// subsequences in order of appearance.
///
/// The input is given in a FASTA format on *standard input* with a single peptide per FASTA header,
/// which may be hardwrapped with newlines. All overlapping k-mers of a peptide are written to
/// *standard output*, separated by newlines. The k-mer length is configurable with the `-k` option,
/// and is 9 by default.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// DAIGDVAKAYKKAG*S
/// $ umgap prot2kmer < input.fa
/// >header1
/// DAIGDVAKA
/// AIGDVAKAY
/// IGDVAKAYK
/// GDVAKAYKK
/// DVAKAYKKA
/// VAKAYKKAG
/// AKAYKKAG*
/// KAYKKAG*S
/// ```
pub struct ProtToKmer {
    /// The k-mer length
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
