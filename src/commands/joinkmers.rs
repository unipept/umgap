//! The `umgap joinkmers` command.

use std::io;
use std::io::Write;
use std::path::PathBuf;

use crate::agg;
use crate::agg::Aggregator;
use crate::errors;
use crate::taxon;
use crate::taxon::TaxonId;
use crate::tree;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Aggregates a TSV stream of peptides and taxon IDs
///
/// The `umgap joinkmers` command takes tab-separated peptides and taxon IDs, aggregates the
/// taxon IDs where consecutive peptides are equal and outputs a tab-separated triple of peptide,
/// consensus taxon ID and taxon rank.
///
/// The input is given on *standard input*. If it is sorted on the first column, a complete mapping
/// from strings to aggregated taxa and its rank will be written to *standard output*. It is
/// meant to be used after an `umgap splitkmers` and `sort`, and it's output is ideal for `umgap
/// buildindex`, but there may be further uses.
///
/// The aggregation strategy used in this command to find a consensus taxon is the hybrid approach
/// of the `umgap taxa2agg` command, with a 95% factor. This keeps the result close to the lowest
/// common ancestor, but filters out some outlying taxa.
///
/// The taxonomy to be used is passed as an argument to this command. This is a preprocessed version
/// of the NCBI taxonomy.
///
/// ```sh
/// $ cat input.tsv
/// AAAAA	34924
/// AAAAA	30423
/// AAAAA	5678
/// BBBBBB	48890
/// BBBBBB	156563
/// $ umgap joinkmers taxons.tsv < input.tsv
/// AAAAA	2759	superkingdom
/// BBBBBB	9153	family
/// ```
pub struct JoinKmers {
    /// An NCBI taxonomy TSV-file as processed by Unipept
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,
}

/// Implements the joinkmers command.
pub fn joinkmers(args: JoinKmers) -> errors::Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let taxons = taxon::read_taxa_file(args.taxon_file)?;

    // Parsing the Taxa file
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let ranksnapping = tree.snapping(&by_id, true);
    let validsnapping = tree.snapping(&by_id, false);
    let aggregator = tree::mix::MixCalculator::new(tree.root, &by_id, 0.95);

    let mut emit = |kmer: &str, tids: Vec<(TaxonId, f32)>| {
        let counts = agg::count(tids.into_iter());
        if let Ok(aggregate) = aggregator.aggregate(&counts) {
            let taxon = ranksnapping[aggregate].unwrap();
            let rank = by_id.get_or_unknown(taxon).unwrap().rank;
            writeln!(handle, "{}\t{}\t{}", kmer, taxon, rank)
        } else {
            Ok(())
        }
    };

    // Iterate over records and emit groups
    let mut current_kmer: Option<String> = Option::None;
    let mut current_tids = vec![];
    for record in reader.deserialize() {
        let (kmer, tid): (String, TaxonId) = record?;
        if let Some(c) = current_kmer {
            if c != kmer {
                emit(&c, current_tids)?;
                current_tids = vec![];
            }
        } else {
            current_tids = vec![];
        }
        current_kmer = Some(kmer);
        if let Some(validancestor) = validsnapping[tid] {
            current_tids.push((validancestor, 1.0));
        }
    }
    if let Some(c) = current_kmer {
        emit(&c, current_tids)?;
    }

    Ok(())
}
