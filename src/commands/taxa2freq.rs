//! The `umgap taxa2freq` command.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::PathBuf;

use crate::errors;
use crate::rank;
use crate::rank::Rank;
use crate::taxon;

#[structopt(verbatim_doc_comment)]
/// Counts ranked taxon occurrences in a stream of taxon IDs
///
/// The `umgap taxa2freq` command creates a frequency table of a list of taxa on a given target rank
/// (species by default).
///
/// The input is given on *standard input*, a single taxon ID on each line. Each taxon that is more
/// specific than the target rank is counted towards its ancestor on the target rank. Each taxon
/// less specific than the target rank is counted towards root. The command outputs a TSV table of
/// counts, taxon IDs and their names.
///
/// The taxonomy to be used is passed as an argument to this command. This is a preprocessed version
/// of the NCBI taxonomy.
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
/// $ umgap taxa2freq taxons.tsv < input.txt
/// 2	1	root
/// 9	9606	Homo sapiens
/// ```
///
/// With the `-r` option, the default species rank can be set to any named rank.
///
/// ```sh
/// $ umgap taxa2freq -r phylum taxons.tsv < input.txt
/// 10	7711	Chordata
/// ```
#[derive(Debug, StructOpt)]
pub struct TaxaToFreq {
    /// The rank to show
    #[structopt(
        short = "r",
        long = "rank",
        default_value = "species",
        possible_values = &Rank::variants()
    )]
    pub rank: Rank,

    /// The minimum frequency to be reported
    #[structopt(short = "f", long = "frequency", default_value = "1")]
    pub min_frequency: usize,

    /// An NCBI taxonomy TSV-file as processed by Unipept
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,

    /// Multiple comparative input files.
    #[structopt(parse(from_os_str))]
    pub input_files: Vec<PathBuf>,
}

/// Implements the taxa2freq command.
pub fn taxa2freq(args: TaxaToFreq) -> errors::Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    if args.rank == rank::Rank::NoRank {
        return Err(errors::ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }
    let numfiles = args.input_files.len();

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping =
        tree.filter_ancestors(|tid| by_id.get(tid).map(|t| t.rank == args.rank).unwrap_or(false));

    // Grab stdout
    let stdout = io::stdout();
    let mut stdout = stdout.lock();

    // Writing headers
    write!(stdout, "taxon id,taxon name")?;
    if numfiles == 0 {
        write!(stdout, ",stdin")?;
    } else {
        for filename in args.input_files.iter() {
            write!(stdout, ",{}", filename.to_string_lossy())?;
        }
    }
    writeln!(stdout)?;

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    if numfiles == 0 {
        let stdin = io::stdin();
        count_file(&snapping, &mut counts, 0, 1, Box::new(stdin.lock()))?;
    } else {
        for (i, file) in args.input_files.iter().enumerate() {
            count_file(
                &snapping,
                &mut counts,
                i,
                numfiles,
                Box::new(BufReader::new(File::open(file)?)),
            )?;
        }
    }

    // Print rows
    for (tid, row) in counts {
        let taxon = by_id
            .get(tid)
            .ok_or("LCA taxon id not in taxon list. Check compatibility with index.")?;
        if row.iter().sum::<usize>() > args.min_frequency {
            write!(stdout, "{},{}", taxon.id, taxon.name)?;
            for count in row {
                write!(stdout, ",{}", count)?;
            }
            writeln!(stdout)?;
        }
    }

    Ok(())
}

fn count_file<T: BufRead>(
    snapping: &Vec<Option<taxon::TaxonId>>,
    counts: &mut HashMap<taxon::TaxonId, Vec<usize>>,
    index: usize,
    numfiles: usize,
    file: T,
) -> errors::Result<()> {
    for line in file.lines() {
        match line?.parse::<taxon::TaxonId>() {
            Ok(taxon) => {
                counts
                    .entry(snapping[taxon].unwrap_or(0))
                    .or_insert(vec![0; numfiles])[index] += 1
            }
            Err(_) => (),
        }
    }
    Ok(())
}
