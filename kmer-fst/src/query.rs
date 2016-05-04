extern crate clap;
extern crate csv;
extern crate fst;

mod errors;
mod fasta;

use std::io;
use std::fs;

use errors::Result;

static K : usize = 7;


fn query(fst_filename: &String, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    let reader = try!(get_reader(query_filename));

    println!("fasta_header,peptide,taxon_id");
    for prot in reader.records() {
        let prot = try!(prot);
        for i in 0..(prot.sequence.len() - K + 1) {
            let kmer = &prot.sequence[i..i + K];
            if let Some(taxon_id) = map.get(kmer) {
                println!("{},{},{}", prot.header, kmer, taxon_id)
            }
        }
    }

    Ok(())
}

fn get_reader(query_filename: &String) -> Result<fasta::Reader<Box<io::Read>>> {
    let reader: Box<io::Read> = match query_filename.as_ref() {
        "-" => Box::new(io::stdin()),
        _   => Box::new(try!(fs::File::open(query_filename)))
    };
    Ok(fasta::Reader::new(reader))
}

fn main() {
    let app = clap::App::new("prot2kmer2lca")
        .arg(clap::Arg::with_name("fst")
             .required(true)
             .help("An FST to query"))
        .arg(clap::Arg::with_name("query file")
             .help("A FASTA formatted file of amino acid sequences. Omit or use '-' to read form stdin"));

    let matches = app.get_matches();

    let fst_filename = String::from(matches.value_of("fst").unwrap());
    let query_filename = String::from(matches.value_of("query file").unwrap_or("-"));

    println!("fst {} query {}", fst_filename, query_filename);

    query(&fst_filename, &query_filename).unwrap();
}
