extern crate unipept;
extern crate clap;
extern crate num_rational;

use std::process;
use std::io;
use std::io::BufRead;

use clap::{Arg, App};
use num_rational::Ratio;

use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::taxon::TaxonId;
use unipept::agg::Aggregator;
use unipept::rmq;
use unipept::tree;

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
                               .conflicts_with("both")
                               .possible_values(&["LCA*", "MRTL", "both", "tree"]))
                      .arg(Arg::with_name("factor")
                               .help("Factor used with aggregate both")
                               .short("f")
                               .long("factor")
                               .takes_value(true))
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
    let aggregation = matches.value_of("aggregate").unwrap_or("LCA*");
    let ranked_only = matches.is_present("ranked");
    let factor      = matches.value_of("factor").unwrap_or("0")
                             .parse::<Ratio<usize>>().unwrap_or_else(|_| {
                                 println!("Error: failed to parse the factor.");
                                 process::exit(2);
                             });

    match aggregation {
        "MTRL" => aggregate(rmq::rtl::RTLCalculator::new(taxons, ranked_only)),
        "LCA*" => aggregate(rmq::lca::LCACalculator::new(taxons, ranked_only)),
        "both" => aggregate(rmq::mix::MixCalculator::new(taxons, ranked_only, factor)),
        "tree" => aggregate(tree::lca::LCACalculator::new(taxons, ranked_only)),
        _      => panic!("Unknown aggregation type.")
    }
}

fn aggregate<T: Aggregator>(aggregator: T) {
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

        let aggregate = aggregator.aggregate(&taxons);
        println!("{},{},{}", aggregate.id, aggregate.name, aggregate.rank);
    }
}
