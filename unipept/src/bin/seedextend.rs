
use std::process;
use std::io;
use std::io::Write;

extern crate clap;
use clap::{Arg, App};

extern crate regex;
use regex::Regex;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::taxon;
use unipept::taxon::TaxonId;
use unipept::io::fasta;

const ABOUT: &'static str = "
Seed and extend.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " seedextend")
        .version(PKG_VERSION)
        .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
        .about(ABOUT)
        .arg(Arg::with_name("min-seed-size")
             .help("The minimum length of equal taxa to count as seed. (default 4)")
             .takes_value(true)
             .short("s")
             .long("min-seed-size"))
        .arg(Arg::with_name("max-gap-size")
             .help("The maximum length of a gap between seeds in an extension. (default 0)")
             .takes_value(true)
             .short("g")
             .long("max-gap-size"))
        .arg(Arg::with_name("taxon-file")
             .help("The taxonomy tsv-file.")
             .index(1)
             .required(true))
        .get_matches();
    main_result(
        matches.value_of("min-seed-size").unwrap_or("4"),
        matches.value_of("max-gap-size").unwrap_or("0"),
        matches.value_of("taxon-file").unwrap(), // required so safe unwrap
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(min_seed_size: &str, max_gap_size: &str, taxon_file: &str) -> Result<(), String> {
    let min_seed_size = try!(min_seed_size.parse::<usize>().map_err(|err| err.to_string()));
    let max_gap_size = try!(max_gap_size.parse::<usize>().map_err(|err| err.to_string()));
    let gap_penalty = 0.5; // TODO: move down
    let taxa = try!(taxon::read_taxa_file(taxon_file).map_err(|err| err.to_string()));

    let tree = taxon::TaxonTree::new(&taxa);
    let by_id = taxon::TaxonList::new_with_unknown(taxa, true);

    let separator = Some(try!(Regex::new(r"\s+").map_err(|err| err.to_string())));
    for record in fasta::Reader::new(io::stdin(), separator, false).records() {
        let record = try!(record.map_err(|err| err.to_string()));
        let taxons = try!(record.sequence.iter()
                                         .map(|s| s.parse::<TaxonId>())
                                         .collect::<Result<Vec<TaxonId>,_>>()
                                         .map_err(|err| err.to_string()));
    }
    Ok(())
}
