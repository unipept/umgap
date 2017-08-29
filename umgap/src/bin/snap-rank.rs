
use std::process;
use std::io;
use std::io::Write;
use std::io::BufRead;
use std::collections::HashMap;

extern crate clap;
use clap::{Arg, App};

extern crate umgap;
use umgap::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};

const ABOUT: &'static str = "
Shows counts of taxa on a specified rank.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " show-tree")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("rank")
                               .help("The rank to show.")
                               .short("r")
                               .long("rank")
                               .takes_value(true))
                      .arg(Arg::with_name("taxon-file")
                               .help("The NCBI taxonomy tsv-file")
                               .index(1)
                               .required(true))
                      .arg(Arg::with_name("with-info")
                               .help("Print additional information about the taxa.")
                               .short("i")
                               .long("with-info"))
                      .get_matches();
    main_result(
        matches.value_of("taxon-file").unwrap(), // required argument, so safe
        matches.value_of("rank").unwrap_or("species"),
        matches.is_present("with-info")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(taxons: &str, rank: &str, with_info: bool) -> Result<(), String> {
    let taxons = try!(taxon::read_taxa_file(taxons).map_err(|err| err.to_string()));
    let rank = try!(rank.parse::<taxon::Rank>().map_err(|err| err.to_string()));
    if rank == taxon::Rank::NoRank {
        return Err("Try snapping to an actual rank.".to_string());
    }

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.filter_ancestors(|tid|
        by_id.get(tid).map(|t| t.rank == rank).unwrap_or(false)
    );

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = try!(try!(line.map_err(|err| err.to_string())).parse::<taxon::TaxonId>().map_err(|err| err.to_string()));
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    // Sorting
    let mut counts = counts.iter().collect::<Vec<_>>();
    counts.sort_by_key(|p| p.1);

    // Printing
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for (tid, count) in counts.into_iter() {
        if with_info {
            let taxon = try!(by_id.get(*tid).ok_or("Missing taxon id".to_string()));
            try!(handle.write_fmt(format_args!("{},{},{},{}\n", tid, taxon.name, taxon.rank, count))
                       .map_err(|err| err.to_string()));
        } else {
            try!(handle.write_fmt(format_args!("{},{}\n", tid, count))
                       .map_err(|err| err.to_string()));
        }
    }

    Ok(())
}

