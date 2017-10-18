
use std::io;
use std::process;
use std::io::Write;

#[macro_use(clap_app, crate_version, crate_authors)]
extern crate clap;

extern crate umgap;
use umgap::dna::Strand;
use umgap::dna::translation::TranslationTable;
use umgap::io::fasta;

fn main() {
    let matches = clap_app!(umgap =>
        (version: crate_version!())
        (author: crate_authors!("\n"))
        (about: "The Unipept Metagenomics Analysis Pipeline")
        (@subcommand translate =>
            (about: "Translates DNA on stdin directly into Amino Acid \
                     Sequences on stdout.")
            (@arg METHIONINE: -m --methionine
                "Replace each start-codon with methionine")
            (@arg ALL_FRAMES: -a --("all-frames") conflicts_with[FRAME]
                "Read and output all 6 frames")
            (@arg FRAME: -f --frame [FRAMES] ...
                "Adds a reading frame (1, 2, 3, 1R, 2R or 3R). If multiple \
                 frames are requested, a bar (|) and the name of the frame \
                 will be appended to the fasta header. (Default: only 1)")
            (@arg TABLE: -t --table [INT]
                "Translation table to use (default: 1)")
            (@arg SHOW_TABLE: -s --("show-table")
                "Print the selected table and exit")
        )
    ).get_matches();

    match matches.subcommand() {
        ("translate", Some(matches)) => translate(
            matches.is_present("methionine"),
            matches.value_of("table").unwrap_or("1"),
            matches.is_present("show-table"),
            matches.values_of("frame").map(Iterator::collect).unwrap_or(
                if matches.is_present("all-frames") { vec!["1", "2", "3", "1R", "2R", "3R"] } else { vec!["1"] }
            )),
        _  => { println!("{}", matches.usage()); Ok(()) }
    }.unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn translate(methionine: bool, table: &str, show_table: bool, frames: Vec<&str>) -> Result<(), String> {
    // Parsing the table
    let table = try!(table.parse::<&TranslationTable>().map_err(|err| err.to_string()));

    // Split on show_tables
    if show_table {
        table.print();
    } else {
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

            let forward = Strand::from(sequence[0].as_bytes());
            let reverse = forward.reversed();
            for &(name, frame, reversed) in &frames {
                let strand = if reversed { &reverse } else { &forward };
                try!(writer.write_record(fasta::Record {
                    header: if !append_name { header.clone() } else { header.clone() + "|" + name },
                    sequence: vec![String::from_utf8(table.translate_frame(methionine, strand.frame(frame))).unwrap()]
                }).map_err(|err| err.to_string()));
            }
        }
    }
    Ok(())
}
