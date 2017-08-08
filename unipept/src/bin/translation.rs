
use std::io;
use std::process;
use std::io::Write;

extern crate clap;
use clap::{App, Arg};

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::translation::{TranslationTable, Codon, invert};
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
                      .arg(Arg::with_name("frame")
                               .help("Adds a reading frame (1, 2, 3, 1R, 2R or 3R). If multiple frames are requested, a bar (|) and the name of the frame will be appended to the fasta header. (Default: only 1)")
                               .short("f")
                               .long("frame")
                               .takes_value(true)
                               .multiple(true))
                      .arg(Arg::with_name("all-frames")
                               .help("Read and output all 6 frames.")
                               .short("a")
                               .long("all-frames")
                               .conflicts_with("frame"))
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
        matches.is_present("show-table"),
        matches.values_of("frame").map(Iterator::collect).unwrap_or(
            if matches.is_present("all-frames") { vec!["1", "2", "3", "1R", "2R", "3R"] } else { vec!["1"] }
        ),
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(methionine: bool, table: &str, show_table: bool, frames: Vec<&str>) -> Result<(), String> {
    // Parsing the table
    let table = try!(table.parse::<&TranslationTable>().map_err(|err| err.to_string()));

    // Split on show_tables
    if show_table {
        table.print(); Ok(())
    } else {
        translate(methionine, table, frames)
    }
}

fn translate(methionine: bool, table: &TranslationTable, frames: Vec<&str>) -> Result<(), String> {
    let translator = |codon| table.translate(methionine, &codon);
    let mut writer = fasta::Writer::new(io::stdout(), "", false);

    // Parsing the frames
    let frames = try!(frames.iter().map(|&frame| match frame {
        "1"  => Ok((frame, 1, false)),
        "2"  => Ok((frame, 2, false)),
        "3"  => Ok((frame, 3, false)),
        "1R" => Ok((frame, 1, true)),
        "2R" => Ok((frame, 2, true)),
        "3R" => Ok((frame, 3, true)),
        _    => return Err(format!("Invalid frame: {}", frame))
    }).collect::<Result<Vec<(&str, usize, bool)>, String>>());
    let append_name = frames.len() > 1;

    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = try!(record.map_err(|err| format!("Something went wrong during the reading: {}", err)));

        let forward = sequence[0].as_bytes().iter().map(|b| *b).collect::<Vec<u8>>();
        let reverse = forward.iter().map(|b| *b).map(invert).rev().collect::<Vec<u8>>();
        for &(name, frame, reversed) in &frames {
            let strand = if reversed { &reverse } else { &forward };
            try!(writer.write_record(fasta::Record {
                header: if !append_name { header.clone() } else { header.clone() + "|" + name },
                sequence: vec![String::from_utf8(strand[frame - 1..].chunks(3).filter(|t| t.len() == 3).map(Codon::from).map(&translator).collect()).unwrap()],
            }).map_err(|err| err.to_string()));
        }
    }

    Ok(())
}
