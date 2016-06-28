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
                      .arg(Arg::with_name("method")
                               .help("The method to use (default RMQ)")
                               .short("m")
                               .long("method")
                               .takes_value(true)
                               .possible_values(&["RMQ", "tree"]))
                      .arg(Arg::with_name("aggregate")
                               .help("The aggregation to use (default LCA*).")
                               .short("a")
                               .long("aggregate")
                               .takes_value(true)
                               .possible_values(&["LCA*", "MRTL", "both"]))
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
    let method      = matches.value_of("method").unwrap_or("RMQ");
    let aggregation = matches.value_of("aggregate").unwrap_or("LCA*");
    let ranked_only = matches.is_present("ranked");
    let factor      = matches.value_of("factor").unwrap_or("0")
                             .parse::<Ratio<usize>>().unwrap_or_else(|err| {
                                 println!("Error: failed to parse the factor.");
                                 println!("Example of a valid factor: 3/4.");
                                 println!("{}", err);
                                 process::exit(2);
                             });

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::taxa_vector_by_id(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    let aggregator: Box<Aggregator> = match (method, aggregation) {
        ("RMQ",  "MTRL") => Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id)),
        ("RMQ",  "LCA*") => Box::new(rmq::lca::LCACalculator::new(tree)),
        ("RMQ",  "both") => Box::new(rmq::mix::MixCalculator::new(tree, factor)),
        ("tree", "LCA*") => Box::new(tree::lca::LCACalculator::new(tree.root, &by_id)),
        ("tree", "both") => Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, factor)),
        _                => {
            println!("Invalid method/aggregation combination: {} and {}", method, aggregation);
            process::exit(3);
        }
    };

    let input = io::BufReader::new(io::stdin());
    for line in input.lines() {
        let line = line.unwrap_or_else(|err| {
            println!("Error: failed to read an input line.");
            println!("{}", err);
            process::exit(4);
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

        let aggregate = aggregator.aggregate(&taxons).unwrap_or_else(|err| {
            println!("Error: failed to aggregate some line: {}", err);
            process::exit(5)
        });
        let taxon = by_id[snapping[aggregate].unwrap()].as_ref().unwrap();
        println!("{},{},{}", taxon.id, taxon.name, taxon.rank);
    }
}

