
use std::process;
use std::io;
use std::io::Write;
use std::io::BufRead;

extern crate clap;
use clap::{Arg, App};

extern crate fst;
use fst::Map;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS, taxon};
use unipept::errors::Result;

const ABOUT: &'static str = "
Looks up each line of input in a given FST index and outputs the result.
Lines starting with a '>' are copied. Lines for which no mapping is found
are ignored.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " pept2lca")
        .version(PKG_VERSION)
        .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
        .about(ABOUT)
        .arg(Arg::with_name("taxon-file")
            .help("The NCBI taxonomy tsv-file")
            .short("t")
            .long("taxa")
            .takes_value(true))
        .arg(Arg::with_name("FST index")
            .help("An FST to query")
            .required(true))
        .get_matches();
    main_result(
        matches.value_of("FST index").unwrap(), // required so safe
        matches.value_of("taxon-file"),
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(fst: &str, taxons: Option<&str>) -> Result<()> {
    let fst = try!(Map::from_path(fst));
    let by_id = taxons.map(|taxons| taxon::read_taxa_file(taxons).map_err(|err| err.to_string()).unwrap())
                      .map(|taxons| taxon::taxa_vector_by_id(taxons));
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let line = try!(line);
        if line.starts_with('>') {
            println!("{}", line);
        } else if let Some(lca) = fst.get(line) {
            if let Some(ref by_id) = by_id {
                let taxon = by_id[lca as usize].as_ref().unwrap();
                println!("{},{},{}", taxon.id, taxon.name, taxon.rank);
            } else {
                println!("{}", lca);
            }
        }
    }
    Ok(())
}
