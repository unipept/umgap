
use std::process;
use std::io;
use std::io::Write;

extern crate clap;
use clap::{Arg, App};

extern crate regex;
use regex::Regex;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::taxon::TaxonId;
use unipept::agg::Aggregator;
use unipept::rmq;
use unipept::tree;
use unipept::io::fasta;
use unipept::agg;
use unipept::errors;

const ABOUT: &'static str = "
Aggregates taxa to a single taxon.
";

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
                               .help("Factor used with aggregate both, in [0.0,1.0]")
                               .short("f")
                               .long("factor")
                               .takes_value(true))
                      .arg(Arg::with_name("separator")
                               .help("Regex to split the taxa (default whitespace)")
                               .short("d")
                               .long("separator")
                               .takes_value(true))
                      .arg(Arg::with_name("lower-bound")
                               .help("The smallest input frequency for a taxon to be included in the aggregation")
                               .short("l")
                               .long("lower-bound")
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
        matches.value_of("separator"),
        matches.is_present("ranked"),
        matches.value_of("factor").unwrap_or("0"),
        matches.is_present("scored"),
        matches.value_of("lower-bound").unwrap_or("0")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}


fn with_score(pair: &String) -> errors::Result<(TaxonId, f32)> {
    let split = pair.split('=').collect::<Vec<_>>();
    if split.len() != 2 { try!(Err("Taxon without score")); }
    Ok((try!(split[0].parse::<TaxonId>()), try!(split[1].parse::<f32>())))
}

fn not_scored(tid: &String) -> errors::Result<(TaxonId, f32)> {
    Ok((try!(tid.parse::<TaxonId>()), 1.0))
}

fn main_result(taxons: &str, method: &str, aggregation: &str, separator: Option<&str>, ranked_only: bool, factor: &str, scored: bool, lower_bound: &str) -> Result<(), String> {
    // Parsing the Taxa file
    let taxons = try!(taxon::read_taxa_file(taxons).map_err(|err| err.to_string()));

    // Parsing the factor
    let factor = try!(factor.parse::<f32>().map_err(|err| err.to_string()));

    // Parsing the factor
    let lower_bound = try!(lower_bound.parse::<f32>().map_err(|err| err.to_string()));

    // Parsing the separator regex
    let separator = try!(Regex::new(separator.unwrap_or(r"\s+"))
                               .map(Some)
                               .map_err(|err| err.to_string()));

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    let aggregator: Result<Box<Aggregator>, String> = match (method, aggregation) {
        ("RMQ",  "MRTL") => Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id))),
        ("RMQ",  "LCA*") => Ok(Box::new(rmq::lca::LCACalculator::new(tree))),
        ("RMQ",  "both") => {
            writeln!(&mut io::stderr(), "Warning: this is a hybrid between LCA/MRTL, not LCA*/MRTL").unwrap();
            Ok(Box::new(rmq::mix::MixCalculator::new(tree, factor)))
        },
        ("tree", "LCA*") => Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id))),
        ("tree", "both") => Ok(Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, factor))),
        _                => Err(format!("Invalid method/aggregation combination: {} and {}", method, aggregation))
    };
    let aggregator = try!(aggregator);

    let parser = if scored { with_score } else { not_scored };

    let mut writer = fasta::Writer::new(io::stdout(), ",", false);

    // Iterate over each read
    for record in fasta::Reader::new(io::stdin(), separator, true).records() {
        // Parse the sequence of LCA's
        let record = try!(record.map_err(|err| err.to_string()));
        let taxons = record.sequence.iter()
                                    .map(parser)
                                    .collect::<Result<Vec<(TaxonId, f32)>,_>>();
        let taxons = try!(taxons.map_err(|err| format!("Error reading taxons ({:?}): {}", record, err)));

        // Create a frequency table of taxons for this read (taking into account the lower bound)
        let counts = agg::count(taxons.into_iter());
        let counts = agg::filter(counts, lower_bound);

        // If we don't have a consensus taxon, leave out the read (i.e. consider undetected)
        if !counts.is_empty() {
            let aggregate = try!(aggregator.aggregate(&counts).map_err(|err| err.to_string()));
            let taxon = by_id.get(snapping[aggregate].unwrap()).unwrap();
            try!(writer.write_record(fasta::Record {
                header: record.header,
                sequence: vec![taxon.id.to_string(), taxon.name.to_string(), taxon.rank.to_string()],
            }).map_err(|err| err.to_string()));
        }
    }
    Ok(())
}

