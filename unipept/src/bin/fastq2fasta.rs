
use std::io;
use std::fs;
use std::process;

extern crate clap;
use clap::{Arg, App};

extern crate umgap;
use umgap::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use umgap::io::fasta;
use umgap::io::fastq;
use umgap::errors::Result;

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

fn open_writer(argument: Option<&str>) -> Result<fasta::Writer<Box<io::Write>>> {
    let output_arg = argument.unwrap_or("-");
    let output: Box<io::Write> = if output_arg == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(try!(fs::File::create(output_arg)))
    };
    Ok(fasta::Writer::new(output, "", false))
}

fn open_reader(argument: &str) -> Result<fastq::Records<Box<io::Read>>> {
    let input: Box<io::Read> = if argument == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(try!(fs::File::open(argument)))
    };
    Ok(fastq::Reader::new(input).records())
}

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " fastq2fasta")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("output")
                               .help("The output file, use '-' or just leave out the option for stdout.")
                               .takes_value(true)
                               .long("output")
                               .short("o"))
                      .arg(Arg::with_name("input")
                               .help("The input files. One of these files may be replaced by - to read from stdin.")
                               .multiple(true)
                               .min_values(1)
                               .index(1))
                      .get_matches();

    let mut writer = open_writer(matches.value_of("output")).unwrap_or_else(|err| {
        println!("Couldn't open writer: {}", err);
        process::exit(1);
    });

    let readers = matches.values_of("input")
                         .ok_or("At least one input file (or -) required.")
                         .unwrap_or_else(|err| {
                             println!("Error when opening the input files: {}", err);
                             process::exit(1)
                         })
                         .map(open_reader)
                         .collect::<Result<Vec<fastq::Records<Box<io::Read>>>>>()
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
                sequence: vec![record.sequence],
            }).unwrap_or_else(|err| {
                println!("Something went wrong during the writing: {}", err);
                process::exit(3)
            });
        }
    }
}
