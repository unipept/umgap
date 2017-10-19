
use std::process;
use std::io;
use std::io::Write;
use std::borrow::Cow;

extern crate clap;
use clap::{Arg, App};

extern crate umgap;
use umgap::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use umgap::io::fasta;
use umgap::errors::Result;

const ABOUT: &'static str = "
Splits each protein sequence in a FASTA format into a list of kmers.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " prot2kmer")
        .version(PKG_VERSION)
        .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
        .about(ABOUT)
        .arg(Arg::with_name("length")
            .help("The k in kmers")
            .short("k")
            .long("length")
            .takes_value(true)
            .default_value("8"))
        .get_matches();
    main_result(matches.value_of("length").unwrap()).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(k: &str) -> Result<()> {
    let k = try!(k.parse::<usize>().map_err(|_| "Couldn't parse kmer length"));

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = try!(record);
        if sequence[0].len() < k { continue }
        try!(writer.write_record(fasta::Record {
            header: header,
            sequence: sequence[0].as_bytes().windows(k)
                                 .map(String::from_utf8_lossy).map(Cow::into_owned)
                                 .collect(),
        }));
    }
    Ok(())
}
