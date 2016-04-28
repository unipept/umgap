
extern crate unipept;

use std::env;
use std::io;
use std::io::BufRead;
use std::fs::File;

use unipept::taxon::Taxon;
use unipept::taxon::TaxonId;
use unipept::lca::LCACalculator;

fn main() {
    //print!("Reading the taxons tree... ");
    let filename = env::args().skip(1).next().unwrap();
    let taxon_file = File::open(filename).unwrap();
    let reader = io::BufReader::new(taxon_file);
    let calculator = LCACalculator::new(reader.lines().map(|mline| {
        let line = mline.unwrap();
        line.parse::<Taxon>().unwrap()
    }).collect());
    //println!("done.");

    //print!("Reading the input taxons... ");
    let input = io::BufReader::new(io::stdin());
    for line in input.lines() {
        let taxons = line.unwrap()
                         .trim_right()
                         .split(' ')
                         .map(|tid| tid.parse::<TaxonId>().unwrap())
                         .collect();
        println!("{}", calculator.calc_lca(&taxons, false));
    }
    //println!("done.");
}

