extern crate csv;
extern crate fst;

mod errors;

use std::fs::File;
use std::io;
use std::io::BufRead;

use errors::Result;


fn query(fst_filename: &String, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    println!("Number of kmers: {}", map.len());

    let query_file = try!(File::open(query_filename));
    let reader = io::BufReader::new(query_file);

    for query_sequence in reader.lines() {
        let query_sequence = try!(query_sequence);
        let taxon_id = map.get(&query_sequence);
        match taxon_id {
            None            => (),  // println!("No match for {}", &query_sequence),
            Some(taxon_idx) => println!("{}", taxon_idx),
        }
    }

    return Ok(());
}

fn main() {
    use std::env;

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        panic!("Please supply an fst and query file with protein sequences.")
    }
    let fst_filename = &args[1];
    let query_filename = &args[2];

    query(fst_filename, query_filename).unwrap();
}
