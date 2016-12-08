
use std::process;
use std::io;
use std::io::Write;

extern crate clap;
use clap::{Arg, App};

extern crate num_rational;
use num_rational::Ratio;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::taxon::TaxonId;
use unipept::agg::Aggregator;
use unipept::rmq;
use unipept::tree;
use unipept::io::fasta;

const ABOUT: &'static str = "
Aggregates taxa to a single taxon.
";

fn main_result(taxons: &str, method: &str, aggregation: &str, ranked_only: bool, factor: &str, scored: bool) -> Result<(), String> {
    // Parsing the Taxa file
    let taxons = try!(taxon::read_taxa_file(taxons).map_err(|err| err.to_string()));

    // Parsing the aggregation method
    let factor = try!(factor.parse::<Ratio<usize>>().map_err(|err| err.to_string()));

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::taxa_vector_by_id(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    let aggregator: Result<Box<Aggregator>, String> = match (method, aggregation) {
        ("RMQ",  "MTRL") => Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id))),
        ("RMQ",  "LCA*") => Ok(Box::new(rmq::lca::LCACalculator::new(tree))),
        ("RMQ",  "both") => Ok(Box::new(rmq::mix::MixCalculator::new(tree, factor))),
        ("tree", "LCA*") => Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id))),
        ("tree", "both") => Ok(Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, factor))),
        _                => Err(format!("Invalid method/aggregation combination: {} and {}", method, aggregation))
    };
    let aggregator = try!(aggregator);

    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = try!(record.map_err(|err| err.to_string()));
        let taxons = record.sequence.trim_right()
                                    .split(' ')
                                    .map(|tid| tid.parse::<TaxonId>())
                                    .collect::<Result<Vec<TaxonId>,_>>();
        let taxons = try!(taxons.map_err(|err| format!("{:?}", record)));
        let aggregate = try!(aggregator.aggregate(&taxons).map_err(|err| err.to_string()));
        let taxon = by_id[snapping[aggregate].unwrap()].as_ref().unwrap();
        println!(">{}", record.header);
        println!("{},{},{}", taxon.id, taxon.name, taxon.rank);
    }
    Ok(())
}

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " taxa2lca")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("scored")
                               .help("Each taxon is followed by a score between 0 and 1")
                               .short("s")
                               .long("scored"))
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
    main_result(
        matches.value_of("taxon-file").unwrap(), // required argument, so safe
        matches.value_of("method").unwrap_or("RMQ"),
        matches.value_of("aggregate").unwrap_or("LCA*"),
        matches.is_present("ranked"),
        matches.value_of("factor").unwrap_or("0"),
        matches.is_present("scored")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

