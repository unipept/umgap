extern crate csv;
extern crate fst;
extern crate itertools;

mod errors;
mod fasta;

use std::str::Chars;

use itertools::PutBackN;

use errors::Result;

static K : usize = 7;


fn query(fst_filename: &String, query_filename: &String) -> Result<()> {
    let map = try!(fst::Map::from_path(fst_filename));
    let reader = try!(fasta::Reader::from_file(query_filename));

    println!("fasta_header,peptide,taxon_id");
    for prot in reader.records() {
        let prot = try!(prot);
        let mut sequence = PutBackN::new(prot.sequence.chars());
        while let Some(kmer) = peek_chars(&mut sequence, K) {
            let taxon_id = map.get(&kmer);
            if let Some(taxon_id) = taxon_id {
                println!("{},{},{}", prot.header, kmer, taxon_id)
            }
            sequence.next();
        }
    }

    Ok(())
}

fn peek_chars(it: &mut PutBackN<Chars>, n: usize) -> Option<String> {
    let mut tmp = String::new();

    for _ in 0..n {
        if let Some(c) = it.next() {
            tmp.push(c);
        } else {
            return None;
        }
    }

    for i in tmp.chars().rev() {
        it.put_back(i);
    }

    Some(tmp)
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
