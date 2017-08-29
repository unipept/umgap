
use std::process;
use std::io;
use std::io::Write;
use std::io::BufRead;
use std::collections::HashMap;
use std::ops::Add;

extern crate clap;
use clap::{Arg, App};

extern crate regex;
use regex::Regex;

#[macro_use]
extern crate serde_json;
use serde_json::value;

extern crate umgap;
use umgap::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use umgap::taxon::TaxonId;
use umgap::agg::Aggregator;
use umgap::rmq;
use umgap::tree::tree;
use umgap::io::fasta;
use umgap::agg;
use umgap::errors;

const ABOUT: &'static str = "
Aggregates taxa to a JSON tree for usage in the unipept-visualisations
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
    main_result(
        matches.value_of("taxon-file").unwrap(), // required argument, so safe
        matches.is_present("ranked")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(taxons: &str, ranked_only: bool) -> Result<(), String> {
    let taxons = try!(taxon::read_taxa_file(taxons).map_err(|err| err.to_string()));

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = try!(try!(line.map_err(|err| err.to_string())).parse::<taxon::TaxonId>().map_err(|err| err.to_string()));
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    let tree = try!(tree::Tree::new(1, &by_id.ancestry(), &counts).map_err(|err| err.to_string()));
    let aggtree = tree.aggregate(&Add::add);
    print!("{}", to_json(&tree, &aggtree, &by_id));

    Ok(())
}

fn to_json(node: &tree::Tree<usize>, aggnode: &tree::Tree<usize>, by_id: &taxon::TaxonList) -> value::Value {
    let root = by_id.get(node.root).unwrap();
    json!({
        "name": root.name,
        "id": node.root,
        "data": {
            "count": aggnode.value,
            "valid_taxon": if root.valid { "1" } else { "0" },
            "rank": format!("{}", root.rank),
            "self_count": node.value
        },
        "children": node.children.iter().zip(aggnode.children.iter())
                                 .map(|(n, s)| to_json(n, s, by_id))
                                 .collect::<Vec<_>>()
    })
}

