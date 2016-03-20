extern crate csv;
extern crate fst;

mod errors;

use std::fs::File;
use std::io;
use std::io::BufRead;

use errors::Result;


fn build_fst(csv_filename: &String, fst_filename: &String) -> Result<()> {
    let csv_file = try!(File::open(csv_filename));
    let mut reader = csv::Reader::from_reader(io::BufReader::new(csv_file))
        .has_headers(false)
        .delimiter(b'\t');

    let writer = io::BufWriter::new(try!(File::create(fst_filename)));
    let mut map = try!(fst::MapBuilder::new(writer));

    for record in reader.decode() {
        let (_id, kmer, lca): (String, String, u64) = try!(record);
        try!(map.insert(kmer, lca));
    }

    try!(map.finish());

    return Ok(());
}

fn query_fst(fst_filename: &String, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    println!("Number of elements: {}", map.len());

    let query_file = try!(File::open(query_filename));
    let reader = io::BufReader::new(query_file);

    for query_sequence in reader.lines() {
        let query_sequence = try!(query_sequence);
        let taxon_id = map.get(&query_sequence);
        match taxon_id {
            None => println!("No match for {}", &query_sequence),
            Some(taxon_idx) => println!("{}", taxon_idx),
        }
    }

    return Ok(());
}

fn main() {
    use std::env;

    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        panic!("Please supply an input, output and query file.")
    }
    let input_filename = &args[1];
    let output_filename = &args[2];
    let query_filename = &args[3];

    build_fst(input_filename, output_filename).unwrap();
    query_fst(output_filename, query_filename).unwrap();
}
