
use std::io;
use std::process;
use std::io::Write;

extern crate clap;
use clap::{App, Arg};

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::translation::{TranslationTable, Codon};
use unipept::io::fasta;

const ABOUT: &'static str = "
Translates DNA on stdin directly into Amino Acid Sequences on stdout.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " translation")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .arg(Arg::with_name("methionine")
                               .help("Replace each start-codon with methionine")
                               .short("m")
                               .long("methionine"))
                      .arg(Arg::with_name("table")
                               .help("Translation table to use (default: 1)")
                               .short("t")
                               .long("table")
                               .takes_value(true))
                      .arg(Arg::with_name("show-table")
                               .help("Print the selected table and exit.")
                               .short("s")
                               .long("show-table"))
                      .get_matches();
    main_result(
        matches.is_present("methionine"),
        matches.value_of("table").unwrap_or("1"),
        matches.is_present("show-table")
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(methionine: bool, table: &str, show_table: bool) -> Result<(), String> {
    // Parsing the table
    let table = try!(table.parse::<&TranslationTable>().map_err(|err| err.to_string()));

    // Split on show_tables
    if show_table {
        table.print(); Ok(())
    } else {
        translate(methionine, table)
    }
}

fn translate(methionine: bool, table: &TranslationTable) -> Result<(), String> {
    let mut writer = fasta::Writer::new(io::stdout(), "", false);
    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = try!(record.map_err(|err| format!("Something went wrong during the reading: {}", err)));
        try!(writer.write_record(fasta::Record {
            header: header,
            sequence: vec![String::from_utf8(sequence[0].as_bytes()
                                 .chunks(3).map(Codon::from)
                                 .map(|c| table.translate(methionine, &c))
                                 .collect()).unwrap()]
        }).map_err(|err| err.to_string()));
    }

    Ok(())
}
