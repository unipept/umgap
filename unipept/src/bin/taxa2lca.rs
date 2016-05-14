
extern crate unipept;
extern crate clap;

use std::process;
use std::io;
use std::io::BufRead;

use clap::{Arg, App};

use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::taxon::TaxonId;
use unipept::lca::LCACalculator;

const ABOUT: &'static str = "
Aggregates taxa to a single taxon.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " taxa2lca")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("ranked")
                               .help("Restrict to taxa with a taxonomic rank")
                               .short("r")
                               .long("ranked"))
                      .arg(Arg::with_name("taxon-file")
                               .help("The NCBI taxonomy tsv-file")
                               .index(1)
                               .required(true))
                      .get_matches();

    let filename = matches.value_of("taxon-file").unwrap(); // required argument
    let taxons   = taxon::read_taxa_file(filename);
    if taxons.is_err() {
        println!("Error: {}", taxons.unwrap_err());
        process::exit(1);
    }
    let calculator = LCACalculator::new(taxons.unwrap());

    let input = io::BufReader::new(io::stdin());
    for line in input.lines() {
        let line = line.unwrap_or_else(|_| {
            println!("Error: failed to read an input line.");
            process::exit(2);
        });

        // Copy FASTA headers
        if line.chars().nth(0) == Some('>') {
            println!("{}", line);
            continue
        }

        let taxons = line.trim_right()
                         .split(' ')
                         .map(|tid| tid.parse::<TaxonId>().unwrap())
                         .collect();
        println!("{}", calculator.calc_lca(&taxons, matches.is_present("ranked")));
    }
}

