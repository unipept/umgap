//! The `umgap splitkmers` command.

use std::io;

use crate::errors;
use crate::taxon::TaxonId;

#[cfg_attr(rustfmt, rustfmt_skip)]
#[structopt(verbatim_doc_comment)]
/// Splits peptides in a CSV stream into K-mers duplicating the records
///
/// The `umgap splitkmers` command takes tab-separated taxon IDs and protein sequences and outputs
/// the k-mers mapped to the taxon IDs.
///
/// The input is given on *standard input* and should be a TSV formatted stream of taxon IDs and a
/// protein sequence from this taxon. The output will be written to *standard output* and consists
/// of a TSV formatted stream of k-mers (*k* is configurable with `-k` and is 9 by default) mapped
/// to the taxa in which they occur.
///
/// This output stream is ready to be grouped by k-mer by sorting and than aggregated into a
/// searchable index, with the `sort`, `umgap joinkmers` and `umgap buildindex` commands.
///
/// ```sh
/// $ cat input.fa
/// 654924  MNAKYDTDQGVGRMLFLGTIGLAVVVGGLMAYGYYYDGKTPSSGTSFHTASPSFSSRYRY
/// 176652  MIKLFCVLAAFISINSACQSSHQQREEFTVATYHSSSICTTYCYSNCVVASQHKGLNVESYTCDKPDPYGRETVCKCTLIKCHDI
/// $ umgap splitkmers < input.fa
/// MNAKYDTDQ       654924
/// NAKYDTDQG       654924
/// AKYDTDQGV       654924
/// KYDTDQGVG       654924
/// YDTDQGVGR       654924
/// ...
/// SPSFSSRYR       654924
/// PSFSSRYRY       654924
/// MIKLFCVLA       176652
/// IKLFCVLAA       176652
/// KLFCVLAAF       176652
/// ...
/// ```
#[derive(Debug, StructOpt)]
pub struct SplitKmers {
    /// The k in k-mers
    #[structopt(short = "k", long = "length", default_value = "9")]
    pub length: usize,

    /// Print only the (k-1)-mer suffixes of the k-mers starting with this character
    #[structopt(short = "p", long = "prefix", default_value = "")]
    pub prefix: String,
}

/// Implements the splitkmers command.
pub fn splitkmers(args: SplitKmers) -> errors::Result<()> {
    let byte = args.prefix.as_bytes().first();

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    for record in reader.deserialize() {
        let (tid, sequence): (TaxonId, String) = record?;
        if sequence.len() < args.length {
            continue;
        }
        for kmer in sequence.as_bytes().windows(args.length) {
            if byte.is_none() {
                writer.serialize((String::from_utf8_lossy(kmer), tid))?;
            } else if *byte.unwrap() == kmer[0] {
                writer.serialize((String::from_utf8_lossy(&kmer[1..]), tid))?;
            }
        }
    }

    Ok(())
}
