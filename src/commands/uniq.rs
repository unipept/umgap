//! The `umgap uniq` command.

use std::io;

use crate::errors;
use crate::io::fasta;

#[structopt(verbatim_doc_comment)]
/// Joins consecutive FASTA records with the same header
///
/// The `umgap uniq` command can for example be used to join together the predictions of 2 paired
/// ends before aggregation.
///
/// The input is given in a FASTA format on *standard input*. The content of all consecutive records
/// with the same FASTA header is joined under a single header, separated by newlines (or another
/// separator set with `-s`).
///
/// A delimiter can be passed with the `-d` option to drop this delimiter and everything after it
/// from the FASTA header before comparing it.
///
/// ```sh
/// $ cat input.fa
/// >header1/1
/// 147206
/// 240495
/// >header1/2
/// 1883
/// 1
/// 1883
/// 1883
/// $ umgap uniq -d / < input.fa
/// >header1
/// 147206
/// 240495
/// 1883
/// 1
/// 1883
/// 1883
/// ```
#[derive(Debug, StructOpt)]
pub struct Uniq {
    /// Separator between output items
    #[structopt(short = "s", long = "separator", default_value = "\n")]
    pub separator: String,

    /// Wrap the output sequences
    #[structopt(short = "w", long = "wrap")]
    pub wrap: bool,

    /// Strip FASTA headers after this string
    #[structopt(short = "d", long = "delimiter")]
    pub delimiter: Option<String>,
}

/// Implements the uniq command.
pub fn uniq(args: Uniq) -> errors::Result<()> {
    let mut last = None::<fasta::Record>;
    let mut writer = fasta::Writer::new(io::stdout(), &args.separator, args.wrap);
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let mut record = record?;
        if let Some(ref delimiter) = args.delimiter {
            record
                .header
                .truncate(record.header.find(delimiter).unwrap_or(record.header.len()));
        }
        if let Some(ref mut rec) = last {
            if rec.header == record.header {
                rec.sequence.extend(record.sequence);
            } else {
                writer.write_record_ref(rec)?;
                *rec = record;
            }
        } else {
            last = Some(record);
        }
    }
    if let Some(rec) = last {
        writer.write_record(rec)?;
    }
    Ok(())
}
