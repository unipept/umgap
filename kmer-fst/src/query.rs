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
    use std::env;

    let mut args: Vec<String> = env::args().collect();

    // If no second argument is given, assume a query is given through stdin
    if args.len() == 2 {
        args.push(String::from("-"))
    }

    if args.len() != 3 {
        panic!("Please either supply an FST and query file with protein \
               sequences in FASTA format, or supply an FST and feed such \
               a FASTA query file through stdin.")
    }

    let fst_filename = &args[1];
    let query_filename = &args[2];

    query(fst_filename, query_filename).unwrap();
}
