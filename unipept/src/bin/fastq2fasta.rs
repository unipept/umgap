
use std::io;
use std::fs;
use std::process;

extern crate clap;
use clap::{Arg, App};

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::io::fasta;
use unipept::io::fastq;
use unipept::errors::Result;

const ABOUT: &'static str = "
Interleaves a number of FASTQ files into a single FASTA file.
";

struct Zip<E, I: Iterator<Item=E>> {
    parts: Vec<I>,
}

impl<E, I: Iterator<Item=E>> Zip<E, I> {
    fn new(parts: Vec<I>) -> Self {
        Zip { parts: parts }
    }
}

impl<E, I: Iterator<Item=E>> Iterator for Zip<E, I> {
    type Item = Vec<E>;
    fn next(&mut self) -> Option<Self::Item> {
        self.parts.iter_mut()
            .map(|part| part.next())
            .collect()
    }
}

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " fastq2fasta")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("output")
                               .help("The output file, use '-' or just leave out the option for stdout.")
                               .takes_value(true)
                               .default_value("-")
                               .long("output")
                               .short("o"))
                      .arg(Arg::with_name("input")
                               .help("The input files. One of these files may be replaced by - to read form stdin.")
                               .multiple(true)
                               .min_values(1)
                               .index(1))
                      .get_matches();

    let mut writer = fasta::Writer::new(matches.value_of("output")
                                               .ok_or(io::Error::new(io::ErrorKind::Other, "oh no"))
                                               .and_then(|filename| fs::File::create(filename))
                                               .unwrap());
    let readers = matches.values_of("input").expect("clap forces force a value")
                         .map(|argument| argument.to_string())
                         .map(|filename| fs::File::open(filename)
                                                  .map_err(From::from)
                                                  .map(fastq::Reader::new)
                                                  .map(fastq::Reader::records))
                         .collect::<Result<Vec<fastq::Records<fs::File>>>>()
                         .unwrap_or_else(|err| {
                             println!("Error when opening the input files: {}", err);
                             process::exit(1)
                         });
    for recordzip in Zip::new(readers) {
        for record in recordzip {
            let record = record.unwrap_or_else(|err| {
                println!("Something went wrong during the reading: {}", err);
                process::exit(2)
            });
            writer.write_record(fasta::Record {
                header: record.header,
                sequence: record.sequence,
            }).unwrap_or_else(|err| {
                println!("Something went wrong during the writing: {}", err);
                process::exit(3)
            });
        }
    }
}
