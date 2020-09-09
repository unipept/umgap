//! The `umgap snaptaxon` command.

use std::io;
use std::io::BufRead;
use std::io::Write;
use std::path::PathBuf;

use crate::errors;
use crate::rank;
use crate::rank::Rank;
use crate::taxon;
use crate::taxon::TaxonId;

#[structopt(verbatim_doc_comment)]
/// Snaps taxon IDs to a specific rank or listed taxa
///
/// The `umgap snaptaxon` command takes one or more taxon IDs. For each taxon, it searches amongst
/// its ancestors if any are of the specified rank (e.g., `-r species`), or are one of the listed
/// taxa (e.g., `-t 1598 -t 1883`). If so, the taxon is replaced by the most specific of these
/// matching taxa. Otherwise, it is mapped to the root of the taxonomy.
///
/// The input is given on *standard input* and may be any sequence of FASTA headers and/or lines
/// containing a single taxon ID. The FASTA headers (if any) are just copied over to *standard
/// output*.
///
/// The taxonomy to be used is passed as an argument to this command. This is a preprocessed version
/// of the NCBI taxonomy.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// 888268
/// 186802
/// 1598
/// 1883
/// $ umgap snaptaxon 2020-04-taxons.tsv -r order < ~/input.fa
/// >header1
/// 38820
/// 186802
/// 186826
/// 85011
/// $ umgap snaptaxon 2020-04-taxons.tsv -t 1239 2 < ~/input.fa
/// >header1
/// 1
/// 1239
/// 1239
/// 2
/// ```
#[derive(Debug, StructOpt)]
pub struct SnapTaxon {
    /// The rank to snap towards.
    #[structopt(short = "r", long = "rank")]
    pub rank: Option<Rank>,

    /// A taxon to snap towards (allow multiple times).
    #[structopt(short = "t", long = "taxons")]
    pub taxons: Vec<TaxonId>,

    /// Include the invalidated taxa from the taxonomy
    #[structopt(short = "i", long = "invalid")]
    pub invalid: bool,

    /// An NCBI taxonomy TSV-file as processed by Unipept
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,
}

/// Implements the snaptaxon command.
pub fn snaptaxon(args: SnapTaxon) -> errors::Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    if args.rank.map(|r| r == rank::Rank::NoRank).unwrap_or(false) {
        return Err(errors::ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = tree.filter_ancestors(|tid| {
        args.taxons.contains(&tid)
            || by_id
                .get(tid)
                .map(|t| {
                    (args.invalid || t.valid) && args.rank.map(|r| t.rank == r).unwrap_or(false)
                })
                .unwrap_or(false)
    });

    // Read and count taxon ranks
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(handle, "{}", line)?;
        } else {
            let taxon = line.parse::<taxon::TaxonId>()?;
            let snapped = snapping[taxon].unwrap_or(0);
            writeln!(handle, "{}", snapped)?;
        }
    }

    Ok(())
}
