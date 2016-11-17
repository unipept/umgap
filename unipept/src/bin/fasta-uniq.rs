
use std::io;
use std::process;

extern crate clap;
use clap::{App, Arg};

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::io::fasta;

const ABOUT: &'static str = "
Concatenates the data strings of all consecutive entries with the same header.
";

fn write_record(writer: &mut fasta::Writer<io::Stdout>, record: &fasta::Record) {
    writer.write_record_ref(record).unwrap_or_else(|err| {
        println!("Error occured while writing records: {}", err);
        process::exit(2);
    });
}

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " fasta-uniq")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("separator")
                               .help("String to separate joined sequences with (default '\\n')")
                               .short("s")
                               .long("separator")
                               .takes_value(true))
                      .arg(Arg::with_name("keep")
                               .help("Keep newlines in the sequence.")
                               .short("k")
                               .long("keep"))
                      .get_matches();
    let separator  = matches.value_of("separator").unwrap_or("\n");
    let keep       = matches.is_present("keep");
    let mut last   = None::<fasta::Record>;
    let mut writer = fasta::Writer::new(io::stdout(), keep);
    for record in fasta::Reader::new(io::stdin(), keep).records() {
        let record = record.unwrap_or_else(|err| {
            println!("Something went wrong during the reading: {}", err);
            process::exit(1);
        });
        if let Some(ref mut rec) = last {
            if rec.header == record.header {
                rec.sequence.extend(separator.chars());
                rec.sequence.extend(record.sequence.chars());
            } else {
                write_record(&mut writer, rec);
                *rec = record;
            }
        } else {
            last = Some(record);
        }
    }
    if let Some(rec) = last {
        write_record(&mut writer, &rec);
    }
}

