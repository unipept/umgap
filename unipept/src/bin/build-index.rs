extern crate csv;
extern crate fst;

use std::fs::File;
use std::io;

extern crate unipept;
use unipept::errors::Result;


fn build(csv_filename: &String, fst_filename: &String) -> Result<()> {
    let csv_file = try!(File::open(csv_filename));
    let mut reader = csv::Reader::from_reader(io::BufReader::new(csv_file))
        .has_headers(false)
        .delimiter(b'\t');

    let writer = io::BufWriter::new(try!(File::create(fst_filename)));
    let mut map = try!(fst::MapBuilder::new(writer));

    for record in reader.decode() {
        let (kmer, lca): (String, u64) = try!(record);
        try!(map.insert(kmer, lca));
    }

    try!(map.finish());

    Ok(())
}

fn main() {
    use std::env;

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        panic!("Please supply an input tsv and output name for the fst.")
    }
    let input_filename = &args[1];
    let output_filename = &args[2];

    build(input_filename, output_filename).unwrap();
}
