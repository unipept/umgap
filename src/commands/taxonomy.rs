//! The `umgap taxonomy` command.

use std::io;
use std::io::BufRead;
use std::io::Write;
use std::path::PathBuf;

use crate::errors;
use crate::rank;
use crate::taxon;

#[cfg_attr(rustfmt, rustfmt_skip)]
#[structopt(verbatim_doc_comment)]
/// Includes info in a stream of taxon IDs
///
/// The `umgap taxonomy` command takes one or more taxon IDs as input, searches for them in a
/// taxonomy and outputs more information about them in a CSV format.
///
/// The input is given on *standard input* and may be any sequence of FASTA headers and/or lines
/// containing a single taxon ID. A CSV header is printed to *standard output*. The FASTA headers
/// (any line starting with a `>`) are just copied over. Each of the taxon IDs on the other lines
/// is looked up in a taxonomy, and the ID, name and rank of the taxon are written out separated by
/// commas.
///
/// A taxonomy file must be passed as argument.
///
/// ```sh
/// $ cat input.fa
/// 2026807
/// 888268
/// 186802
/// 1598
/// 1883
/// $ umgap taxonomy taxons.tsv < input.fa
/// taxon_id,taxon_name,taxon_rank
/// 2026807,Zetaproteobacteria bacterium,species
/// 888268,Dichanthelium oligosanthes,species
/// 186802,Clostridiales,order
/// 1598,Lactobacillus reuteri,species
/// 1883,Streptomyces,genus
/// ```
///
/// The `-H` flag can be used to suppress the CSV header, for instance when dealing with FASTA input.
///
/// ```sh
/// $ cat input2.fa
/// >header1
/// 2026807
/// 888268
/// 186802
/// 1598
/// 1883
/// $ umgap taxonomy -H taxons.tsv < input2.fa
/// >header1
/// 2026807,Zetaproteobacteria bacterium,species
/// 888268,Dichanthelium oligosanthes,species
/// 186802,Clostridiales,order
/// 1598,Lactobacillus reuteri,species
/// 1883,Streptomyces,genus
/// ```
///     
/// The `-a` flag can be used to request a complete ranked lineage.
///
/// ```sh
/// $ cat input3.fa
/// 888268
/// $ umgap taxonomy -a taxons.tsv < input3.fa
/// taxon_id,taxon_name,taxon_rank,superkingdom_id,superkingdom_name,kingdom_id,kingdom_name,subkingdom_id,subkingdom_name,superphylum_id,superphylum_name,phylum_id,phylum_name,subphylum_id,subphylum_name,superclass_id,superclass_name,class_id,class_name,subclass_id,subclass_name,infraclass_id,infraclass_name,superorder_id,superorder_name,order_id,order_name,suborder_id,suborder_name,infraorder_id,infraorder_name,parvorder_id,parvorder_name,superfamily_id,superfamily_name,family_id,family_name,subfamily_id,subfamily_name,tribe_id,tribe_name,subtribe_id,subtribe_name,genus_id,genus_name,subgenus_id,subgenus_name,species_group_id,species_group_name,species_subgroup_id,species_subgroup_name,species_id,species_name,subspecies_id,subspecies_name,varietas_id,varietas_name,forma_id,forma_name
/// 888268,Dichanthelium oligosanthes,species,2759,Eukaryota,33090,Viridiplantae,,,,,35493,Streptophyta,131221,Streptophytina,,,3398,Magnoliopsida,1437197,Petrosaviidae,,,,,38820,Poales,,,,,,,,,4479,Poaceae,147369,Panicoideae,147428,Paniceae,1648011,Dichantheliinae,161620,Dichanthelium,,,,,,,888268,Dichanthelium oligosanthes,,,,,,
/// ```
#[derive(Debug, StructOpt)]
pub struct Taxonomy {
    /// An NCBI taxonomy TSV-file as processed by Unipept
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,

    /// Show the full lineage of a taxon. Ranks below the given taxon
    /// whill be empty.
    #[structopt(short = "a", long = "all")]
    pub all_ranks: bool,

    /// Do not output the CSV header
    #[structopt(short = "H", long = "no-header")]
    pub no_header: bool,
}

/// Implements the taxonomy command.
pub fn taxonomy(args: Taxonomy) -> errors::Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    let by_id = taxon::TaxonList::new(taxons);

    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if !args.no_header {
        write!(handle, "taxon_id,taxon_name,taxon_rank")?;
        if args.all_ranks {
            for rank in rank::Rank::ranks() {
                let rank_name = rank.to_string().replace(" ", "_");
                write!(handle, ",{}_id,{}_name", rank_name, rank_name)?;
            }
        }
        write!(handle, "\n")?;
    }

    for line in stdin.lock().lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(handle, "{}", line)?;
        } else {
            let id = line.parse::<taxon::TaxonId>()?;
            // Map to root if id not found
            let taxon = by_id.get_or_unknown(id)?;
            write!(handle, "{},{},{}", taxon.id, taxon.name, taxon.rank)?;
            if args.all_ranks {
                let lineage = by_id.lineage(id)?;
                for rank in rank::Rank::ranks() {
                    if let Some(l_taxon) = &lineage[rank] {
                        write!(handle, ",{},{}", l_taxon.id, l_taxon.name)?;
                    } else {
                        write!(handle, ",,")?;
                    }
                }
            }
            write!(handle, "\n")?;
        }
    }
    Ok(())
}
