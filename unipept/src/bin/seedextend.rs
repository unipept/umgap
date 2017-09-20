
use std::process;
use std::io;
use std::io::Write;
use std::cmp;

extern crate clap;
use clap::{Arg, App};

extern crate regex;
use regex::Regex;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::taxon;
use unipept::taxon::TaxonId;
use unipept::io::fasta;

const ABOUT: &'static str = "
Seed and extend.
";

fn main() {
    let matches = App::new(PKG_NAME.to_string() + " seedextend")
        .version(PKG_VERSION)
        .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
        .about(ABOUT)
        .arg(Arg::with_name("min-seed-size")
             .help("The minimum length of equal taxa to count as seed. (default 4)")
             .takes_value(true)
             .short("s")
             .long("min-seed-size"))
        .arg(Arg::with_name("max-gap-size")
             .help("The maximum length of a gap between seeds in an extension. (default 0)")
             .takes_value(true)
             .short("g")
             .long("max-gap-size"))
        .arg(Arg::with_name("taxon-file")
             .help("The taxonomy tsv-file.")
             .index(1)
             .required(true))
        .get_matches();
    main_result(
        matches.value_of("min-seed-size").unwrap_or("4"),
        matches.value_of("max-gap-size").unwrap_or("0"),
        matches.value_of("taxon-file").unwrap(), // required so safe unwrap
    ).unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result(min_seed_size: &str, max_gap_size: &str, taxon_file: &str) -> Result<(), String> {
    let min_seed_size = try!(min_seed_size.parse::<usize>().map_err(|err| err.to_string()));
    let max_gap_size = try!(max_gap_size.parse::<usize>().map_err(|err| err.to_string()));
    let _gap_penalty = 0.5; // TODO: move down
    let _taxa = try!(taxon::read_taxa_file(taxon_file).map_err(|err| err.to_string())); // TODO use it

    let mut writer = fasta::Writer::new(io::stdout(), ",", false);

    let separator = Some(try!(Regex::new(r"\s+").map_err(|err| err.to_string())));
    for record in fasta::Reader::new(io::stdin(), separator, false).records() {
        let record = try!(record.map_err(|err| err.to_string()));
        let taxons = try!(record.sequence.iter()
                                         .map(|s| s.parse::<TaxonId>())
                                         .collect::<Result<Vec<TaxonId>,_>>()
                                         .map_err(|err| err.to_string()));
        if taxons.len() == 0 {
            return Err("Frame should contain at least one taxon".to_string());
        }

        let mut seeds = Vec::new();
        let mut start = 0;
        let mut end = 0;
        let mut last_tid = 0;
        let mut same_tid = 0;
        let mut same_max = 0;
        while end < taxons.len() {
            // same tid as last, add to seed.
            if last_tid == taxons[end] {
                same_tid += 1;
                end += 1;
                continue;
            }

            // our gap just became to big
            if last_tid == 0 && same_tid > max_gap_size {
                // add extended seed
                if same_max >= min_seed_size { seeds.push((start, end - same_tid)) }
                start = end;
                last_tid = taxons[end];
                same_tid = 1;
                same_max = 1;
                end += 1;
                continue;
            }

            // don't start with a missing taxon
            if last_tid == 0 && (end - start) == same_tid {
                end += 1;
                start = end;
                continue;
            }

            // another taxon
            if last_tid != 0 { same_max = cmp::max(same_max, same_tid); }
            last_tid = taxons[end];
            same_tid = 1;
            end += 1;
        }
        if same_max >= min_seed_size {
            if last_tid == 0 { end -= same_tid }
            seeds.push((start, end))
        }

        // write it
        try!(writer.write_record(fasta::Record {
            header: record.header,
            sequence: seeds.into_iter()
                           .flat_map(|(start, end)| taxons[start..end].into_iter())
                           .map(|t| t.to_string())
                           .collect()
        }).map_err(|err| err.to_string()));
    }
    Ok(())
}
