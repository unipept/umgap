
use std::io;
use std::process;
use std::io::Write;

extern crate clap;
use clap::{App, Arg};

extern crate regex;
use regex::Regex;

extern crate umgap;
use umgap::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use umgap::io::fasta;

const ABOUT: &'static str = "
Concatenates the data strings of all consecutive entries with the same header.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " fasta-uniq")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("separator")
                               .help("Separator between output items (default the empty string)")
                               .short("s")
                               .long("separator")
                               .takes_value(true))
                      .arg(Arg::with_name("input-separator")
                               .help("Separator regex between input items (default same as separator)")
                               .short("i")
                               .long("input-separator")
                               .takes_value(true))
                      .arg(Arg::with_name("keep")
                               .help("Keep newlines in the input sequence.")
                               .short("k")
                               .long("keep"))
                      .arg(Arg::with_name("wrap")
                               .help("Wrap the output sequences.")
                               .short("w")
                               .long("wrap"))
                      .get_matches();
    main_result(
        matches.value_of("separator").unwrap_or(""),
        matches.value_of("input-separator"),
        matches.is_present("keep"),
        matches.is_present("wrap")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(separator: &str, input_separator: Option<&str>, keep: bool, wrap: bool) -> Result<(), String> {
    // Parsing the input separator regex
    let input_separator = try!(Regex::new(input_separator.unwrap_or(&regex::escape(separator)))
                                     .map_err(|err| err.to_string()));

    let mut last   = None::<fasta::Record>;
    let mut writer = fasta::Writer::new(io::stdout(), separator, wrap);
    for record in fasta::Reader::new(io::stdin(), Some(input_separator), !keep).records() {
        let record = try!(record.map_err(|err| format!("Something went wrong during the reading: {}", err)));
        if let Some(ref mut rec) = last {
            if rec.header == record.header {
                rec.sequence.extend(record.sequence);
            } else {
                try!(writer.write_record_ref(rec)
                           .map_err(|err| format!("Error occured while writing records: {}", err)));
                *rec = record;
            }
        } else {
            last = Some(record);
        }
    }
    if let Some(rec) = last {
        try!(writer.write_record(rec)
                   .map_err(|err| format!("Error while writing last record: {}", err)));
    }
    Ok(())
}

