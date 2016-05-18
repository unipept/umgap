extern crate unipept;
extern crate clap;

use std::process;
use std::io;
use std::io::BufRead;

use clap::{Arg, App};

use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::taxon::TaxonId;
use unipept::agg::Aggregator;
use unipept::lca::LCACalculator;
use unipept::rtl::RTLCalculator;

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
                      .arg(Arg::with_name("aggregate")
                               .help("The aggregation to use (default LCA*).")
                               .short("a")
                               .long("aggregate")
                               .takes_value(true)
                               .conflicts_with("mrtl")
                               .possible_values(&["LCA*", "MRTL"]))
                      .arg(Arg::with_name("mrtl")
                               .help("Short for --aggregate=MRTL")
                               .short("m")
                               .long("mrtl")
                               .conflicts_with("aggregate"))
                      .arg(Arg::with_name("taxon-file")
                               .help("The NCBI taxonomy tsv-file")
                               .index(1)
                               .required(true))
                      .get_matches();

    // Parsing the Taxa file
    let filename = matches.value_of("taxon-file").unwrap(); // required argument
    let taxons = taxon::read_taxa_file(filename).unwrap_or_else(|err| {
        println!("Error: {}", err);
        process::exit(1);
    });

    // Parsing the aggregation method
    let aggregation = if matches.is_present("mrtl") {
        "MTRL"
    } else {
        matches.value_of("aggregate").unwrap_or("LCA*")
    };

    let ranked_only = matches.is_present("ranked");

    match aggregation {
        "MTRL" => aggregate(RTLCalculator::new(taxons), ranked_only),
        "LCA*" => aggregate(LCACalculator::new(taxons), ranked_only),
        _      => panic!("Unknown aggregation type.")
    }
}

fn aggregate<T: Aggregator>(aggregator: T, ranked_only: bool) {
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

        let aggregate = aggregator.aggregate(&taxons, ranked_only);
        println!("{},{},{}", aggregate.id, aggregate.name, aggregate.rank);
    }
}
