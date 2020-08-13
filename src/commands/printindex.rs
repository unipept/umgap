//! The `umgap printindex` command.

use std::io;
use std::path::PathBuf;

use fst::Streamer;

use crate::errors;

#[structopt(verbatim_doc_comment)]
/// Print the key/value pairs in an FST index in TSV format, mostly for debugging.
///
/// ```sh
/// $ umgap printindex tryptic.index
/// ...
/// AAAAADRPANEIGGR	293089
/// AAAAADRPAPAGHDHQAVAR	156981
/// AAAAADRPASQIVR	536018
/// AAAAADRPE	1707
/// AAAAADRPEVHALALR	1883427
/// AAAAADRPFVAEPAR	41275
/// AAAAADRPIAAHAEDESLVR	33010
/// AAAAADRPIR	1988
/// AAAAADRPLAEHGGPVPR	1827
/// ...
/// ```
#[derive(Debug, StructOpt)]
pub struct PrintIndex {
    /// An FST to query
    #[structopt(parse(from_os_str))]
    pub fst_file: PathBuf,
}

/// Implements the printindex command.
pub fn printindex(args: PrintIndex) -> errors::Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    let index = unsafe { fst::Map::from_path(args.fst_file) }?;
    let mut stream = index.stream();

    while let Some((k, v)) = stream.next() {
        writer.serialize((String::from_utf8_lossy(k), v))?;
    }

    Ok(())
}
