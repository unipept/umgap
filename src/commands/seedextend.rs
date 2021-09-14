//! The `umgap seedextend` command.

use std::cmp;
use std::io;
use std::path::PathBuf;

use crate::errors;
use crate::io::fasta;
use crate::taxon;
use crate::taxon::TaxonId;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Selects promising regions in sequences of taxon IDs
///
/// The `umgap seedextend` command takes one or more sequences of taxon IDs and selects regions of
/// consecutive predictions. It can be used to filter out accidental matches of incorrect taxa.
///
/// The input is given in a FASTA format on *standard input*. It should consist of taxon IDs
/// separated by newlines, and the order of these taxa should reflect their location on a peptide,
/// such as output by the `umgap prot2kmer2lca -o` command. As such, 3 consecutive equal IDs
/// representing 9-mers, for instance, indicate a 11-mer match. This so-called seed could still be
/// extended with other taxa, forming an extended seed. The command writes all taxa in any of these
/// extended seeds to *standard output*.
///
/// ```sh
/// $ cat dna.fa
/// >header1
/// CGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCC
/// $ umgap translate -n -a < dna.fa | umgap prot2kmer2lca 9mer.index > input.fa
/// >header1|1
/// 9606 9606 2759 9606 9606 9606 9606 9606 9606 9606 8287
/// >header1|2
/// 2026807 888268 186802 1598 1883
/// >header1|3
/// 1883
/// >header1|1R
/// 27342 2759 155619 1133106 38033 2
/// >header1|2R
/// >header1|3R
/// 2951
/// $ umgap seedextend < input.fa
/// >header1|1
/// 9606 9606 2759 9606 9606 9606 9606 9606 9606 9606 8287
/// >header1|2
/// >header1|3
/// >header1|1R
/// >header1|2R
/// >header1|3R
/// ```
///
/// Taxon IDs are separated by newlines in the actual output, but are separated by spaces in this
/// example.
///
/// The number of consecutive equal IDs to start a seed is 2 by default, and can be changed using
/// the `-s` option. The maximum length of gaps between seeds to join in an extension can be set
/// with `-g`, no gaps are allowed by default.
///
/// The command can be altered to print only the extended seed with the highest score among all
/// extended seeds. Pass a taxonomy using the `-r taxon.tsv` option to activate this. In this scored
/// mode, extended seeds with gaps are given a penalty of 5, which can be made more or less severe
/// (higher or lower) with the `-p` option.
pub struct SeedExtend {
    /// The minimum length of equal taxa to count as seed
    #[structopt(short = "s", long = "min-seed-size", default_value = "2")]
    pub min_seed_size: usize,

    /// The maximum length of a gap between seeds in an extension
    #[structopt(short = "g", long = "max-gap-size", default_value = "0")]
    pub max_gap_size: usize,

    /// Use taxon ranks in given NCBI taxonomy tsv-file to pick extended seed with highest score
    #[structopt(short = "r", long = "ranked", parse(from_os_str))]
    pub ranked: Option<PathBuf>,

    /// The score penalty for gaps in extended seeds
    #[structopt(short = "p", long = "penalty", default_value = "5")]
    pub penalty: usize,
}

/// Implements the seedextend command.
pub fn seedextend(args: SeedExtend) -> errors::Result<()> {
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    let by_id = if let Some(ref tf) = args.ranked {
        let taxa = taxon::read_taxa_file(tf)?;
        Some(taxon::TaxonList::new_with_unknown(taxa, true))
    } else {
        None
    };

    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        let mut taxons = record
            .sequence
            .iter()
            .map(|s| s.parse::<TaxonId>())
            .collect::<std::result::Result<Vec<TaxonId>, _>>()?;
        taxons.push(0);

        let mut seeds = Vec::new();
        let mut start = 0;
        let mut end = 1;
        let mut last_tid = taxons[start];
        let mut same_tid = 1;
        let mut same_max = 1;
        while end < taxons.len() {
            // same tid as last, add to seed.
            if last_tid == taxons[end] {
                same_tid += 1;
                end += 1;
                continue;
            }

            // our gap just became to big
            if last_tid == 0 && same_tid > args.max_gap_size {
                // add extended seed
                if same_max >= args.min_seed_size {
                    seeds.push((start, end - same_tid))
                }
                start = end;
                last_tid = taxons[end];
                same_tid = 1;
                same_max = 1;
                end += 1;
                continue;
            }

            // don't start with a missing taxon
            if last_tid == 0 && (end - start) == same_tid {
                end += 1;
                start = end;
                continue;
            }

            // another taxon
            if last_tid != 0 {
                same_max = cmp::max(same_max, same_tid);
            }
            last_tid = taxons[end];
            same_tid = 1;
            end += 1;
        }
        if same_max >= args.min_seed_size {
            if last_tid == 0 {
                end -= same_tid
            }
            seeds.push((start, end))
        }

        if let Some(ref by_id) = by_id {
            seeds = seeds
                .into_iter()
                .max_by_key(|(s, e)| {
                    taxons
                        .iter()
                        .skip(*s)
                        .take(e - s)
                        .map(|t| by_id.score(*t).unwrap_or(args.penalty))
                        .sum::<usize>()
                })
                .into_iter()
                .collect::<Vec<(usize, usize)>>();
        }

        // write it
        writer
            .write_record(fasta::Record {
                header: record.header,
                sequence: seeds
                    .into_iter()
                    .flat_map(|(start, end)| taxons[start..end].iter())
                    .map(|t| t.to_string())
                    .collect(),
            })
            .map_err(|err| err.to_string())?;
    }
    Ok(())
}
