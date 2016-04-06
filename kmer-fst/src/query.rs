extern crate csv;
extern crate fst;

mod errors;
mod fasta;

use errors::Result;


fn query(fst_filename: &String, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    let reader = try!(fasta::Reader::from_file(query_filename));

    for query in reader.records() {
        let query = try!(query);
        print!("{}", query.header);
        let taxon_id = map.get(&query.sequence[..7]);
        match taxon_id {
            None            => println!("0"),
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
