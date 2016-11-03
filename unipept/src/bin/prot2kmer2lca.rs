#[macro_use] extern crate clap;
extern crate fst;
extern crate itertools;

use std::io;
use std::io::Write;
use std::fs;

extern crate unipept;
use unipept::errors::Error;
use unipept::errors::Result;
use unipept::io::fasta;

/// Reads all the records in a specified FASTA file and queries each in an FST for the LCA's.
/// The result is printed to standard output.
///
/// # Arguments
/// * `fst_filename`: the file which stores the FST
/// * `k`: specifies the length of the k-mers
/// * `query_filename`: the FASTA file containing the records to get the LCA from.
fn query(fst_filename: &String, k: usize, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    let reader = try!(get_reader(query_filename));

    for prot in reader.records() {
        let prot = try!(prot);

        if prot.sequence.len() < k {
            continue
        }

        let lcas = (0..(prot.sequence.len() - k + 1))
            .map(|i| &prot.sequence[i..i + k])
            .filter_map(|kmer| map.get(kmer))
            .map(|lca| lca.to_string())
            .collect::<Vec<_>>()
            .join(" ");

        if ! lcas.is_empty() {
            if let Err(e) = writeln!(io::stdout(), "{}\n{}", prot.header, lcas) {
                if e.kind() == io::ErrorKind::BrokenPipe {
                    break
                } else {
                    return Err(Error::Io(e))
                }
            }
        }
    }

    Ok(())
}

/// Returns the Reader for the given filename, or for stdin if the filename is "-".
///
/// # Arguments
/// * `query_filename`: the file which stores the FST
fn get_reader(query_filename: &String) -> Result<fasta::Reader<Box<io::Read>>> {
    let reader: Box<io::Read> = match query_filename.as_ref() {
        "-" => Box::new(io::stdin()),
        _   => Box::new(try!(fs::File::open(query_filename)))
    };
    Ok(fasta::Reader::new(reader, false))
}

fn main() {
    let app = clap::App::new("prot2kmer2lca")
        .arg(clap::Arg::with_name("fst")
             .required(true)
             .help("An FST to query"))
        .arg(clap::Arg::with_name("k-mer length")
             .required(true)
             .help("The length of the k-mers in the FST"))
        .arg(clap::Arg::with_name("query file")
             .help("A FASTA formatted file of amino acid sequences. \
                   Omit or use '-' to read form stdin"));

    let matches = app.get_matches();

    let fst_filename = String::from(matches.value_of("fst").unwrap());
    let k = value_t!(matches, "k-mer length", usize).unwrap();
    let query_filename = String::from(matches.value_of("query file").unwrap_or("-"));

    query(&fst_filename, k, &query_filename).unwrap();
}
