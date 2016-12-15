
use std::process;
use std::io;
use std::io::Write;
use std::io::BufRead;

extern crate clap;
use clap::{Arg, App};

extern crate fst;
use fst::Map;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
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
        .arg(Arg::with_name("FST index")
            .help("An FST to query")
            .required(true))
        .get_matches();
    main_result(matches.value_of("FST index").unwrap() /* required so safe */)
        .unwrap_or_else(|err| {
            writeln!(&mut io::stderr(), "{}", err).unwrap();
            process::exit(1);
        });
}

fn main_result(fst: &str) -> Result<()> {
    let fst = try!(Map::from_path(fst));
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let line = try!(line);
        if line.starts_with('>') {
            println!("{}", line);
        } else if let Some(lca) = fst.get(line) {
            println!("{}", lca);
        }
    }
    Ok(())
}
