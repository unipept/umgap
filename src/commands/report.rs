//! The `umgap report` command.

use std::collections::HashMap;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::path::PathBuf;

use crate::errors;
use crate::rank;
use crate::rank::Rank;
use crate::taxon;

#[structopt(verbatim_doc_comment)]
/// The `umgap report` command creates a frequency table of a list of taxa.
///
/// The input is given on *standard input*, a single taxon ID on each line. Each taxon at least as
/// specific as a species is counted towards that species, each taxon above species rank is counted
/// as root. The command outputs a CSV table of counts, species taxon IDs and their names.
///
/// For operation, the command requires to be passed a small database created from the uniprot
/// taxonomy as argument.
///
/// ```sh
/// $ cat input.txt
/// 9606
/// 9606
/// 2759
/// 9606
/// 9606
/// 9606
/// 9606
/// 9606
/// 9606
/// 9606
/// 8287
/// $ umgap report taxons.tsv < input.txt
/// 2,1,root
/// 9,9606,Homo sapiens
/// ```
///
/// With the `-r` option, the default species rank can be set to any named rank.
///
/// ```sh
/// $ umgap report -r phylum taxons.tsv < input.txt
/// 10,7711,Chordata
/// ```
#[derive(Debug, StructOpt)]
pub struct Report {
    /// The rank to show
    #[structopt(short = "r", long = "rank", default_value = "species")]
    pub rank: Rank,

    /// The minimum frequency to be reported
    #[structopt(short = "f", long = "frequency", default_value = "1")]
    pub min_frequency: usize,

    /// The NCBI taxonomy tsv-file
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,
}

/// Implements the report command.
pub fn report(args: Report) -> errors::Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    if args.rank == rank::Rank::NoRank {
        return Err(errors::ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping =
        tree.filter_ancestors(|tid| by_id.get(tid).map(|t| t.rank == args.rank).unwrap_or(false));

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = line?.parse::<taxon::TaxonId>()?;
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    let mut counts = counts
        .into_iter()
        .filter(|&(_taxon, count)| count > args.min_frequency)
        .map(|(taxon, count)| (count, taxon))
        .collect::<Vec<(taxon::TaxonId, usize)>>();
    counts.sort();

    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    for (count, taxon) in counts {
        let taxon = by_id
            .get(taxon)
            .ok_or("LCA taxon id not in taxon list. Check compatibility with index.")?;
        writeln!(stdout, "{},{},{}", count, taxon.id, taxon.name)?;
    }

    Ok(())
}
