
use std::io;
use std::process;
use std::io::Write;
use std::io::BufRead;

#[macro_use(clap_app, crate_version, crate_authors)]
extern crate clap;

extern crate fst;

extern crate regex;

extern crate umgap;
use umgap::dna::Strand;
use umgap::dna::translation::TranslationTable;
use umgap::io::fasta;
use umgap::taxon;
use umgap::errors::{Result, ErrorKind};
use umgap::taxon::TaxonId;
use umgap::agg;
use umgap::rmq;
use umgap::tree;

fn main() {
    let matches = clap_app!(umgap =>
        (version: crate_version!())
        (author: crate_authors!("\n"))
        (about: "The Unipept Metagenomics Analysis Pipeline")
        (@subcommand translate =>
            (about: "Translates DNA on stdin directly into Amino Acid \
                     Sequences on stdout.")
            (@arg methionine: -m --methionine
                "Replace each start-codon with methionine")
            (@arg all_frames: -a --("all-frames") conflicts_with[FRAME]
                "Read and output all 6 frames")
            (@arg frame: -f --frame [FRAMES] ...
                "Adds a reading frame (1, 2, 3, 1R, 2R or 3R). If multiple \
                 frames are requested, a bar (|) and the name of the frame \
                 will be appended to the fasta header. (Default: only 1)")
            (@arg table: -t --table [INT]
                "Translation table to use (default: 1)")
            (@arg show_table: -s --("show-table")
                "Print the selected table and exit")
        )
        (@subcommand pept2lca =>
            (about: "Looks up each line of input in a given FST index and \
                     outputs the result. Lines starting with a '>' are \
                     copied. Lines for which no mapping is found are ignored. \
                     Prints more information if the taxa are supplied.")
            (@arg taxon_file: -t --taxa [FILE]
                "The NCBI taxonomy tsv-file")
            (@arg with_input: -i --("with-input")
                "Print the identified peptide along with the output")
            (@arg fst_index: +required "An FST to query")
        )
        (@subcommand prot2kmer2lca =>
            (about: "Reads all the records in a specified FASTA file and \
                     queries the k-mers in an FST for the LCA's.")
            (@arg length: -k --length +required
                "The length of the k-mers in the FST")
            (@arg fst_index: +required "An FST to query")
        )
        (@subcommand taxa2agg =>
            (about: "Aggregates taxa to a single taxon.")
            (@arg scored: -s --scored
                "Each taxon is followed by a score between 0 and 1")
            (@arg ranked: -r --ranked
                "Restrict to taxa with a taxonomic rank")
            (@arg method: -m --method [RMQtree]
                "The method to use for aggregation")
            (@arg aggregate: -a --aggregate [LCAMRTLhybrid]
                "The aggregation to use")
            (@arg factor: -f --factor [RATIO]
                "The factor for the hybrid aggregation, from 0.0 (MRTL) to \
                 1.0 (LCA*)")
            (@arg delimiter: -d --delimiter [REGEX]
                "Regex to split input taxa, default is whitespace")
            (@arg lower_bound: -l --("lower-bound") [INT]
                "The smallest input frequency for a taxon to be included in \
                 the aggregation")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
    ).get_matches();

    match matches.subcommand() {
        ("translate", Some(matches)) => translate(
            matches.is_present("methionine"),
            matches.value_of("table").unwrap_or("1"),
            matches.is_present("show_table"),
            matches.values_of("frame").map(Iterator::collect).unwrap_or(
                if matches.is_present("all_frames") { vec!["1", "2", "3", "1R", "2R", "3R"] } else { vec!["1"] }
            )),
        ("pept2lca", Some(matches)) => pept2lca(
            matches.value_of("fst_index").unwrap(), // required so safe
            matches.value_of("taxon_file"),
            matches.is_present("with_input")),
        ("prot2kmer2lca", Some(matches)) => prot2kmer2lca(
            matches.value_of("fst_index").unwrap(), // required so safe
            matches.value_of("length").unwrap_or("9")),
        ("taxa2agg", Some(matches)) => taxa2agg(
            matches.value_of("taxon_file").unwrap(), // required argument, so safe
            matches.value_of("method").unwrap_or("RMQ"),
            matches.value_of("aggregate").unwrap_or("LCA*"),
            matches.value_of("delimiter"),
            matches.is_present("ranked"),
            matches.value_of("factor").unwrap_or("0"),
            matches.is_present("scored"),
            matches.value_of("lower_bound").unwrap_or("0")),
        _  => { println!("{}", matches.usage()); Ok(()) }
    }.unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn translate(methionine: bool, table: &str, show_table: bool, frames: Vec<&str>) -> Result<()> {
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
            _    => Err(ErrorKind::InvalidInvocation(format!("{} is not a frame", frame)).into())
        }).collect::<Result<Vec<(&str, usize, bool)>>>());
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

fn pept2lca(fst: &str, taxons: Option<&str>, with_input: bool) -> Result<()> {
    let fst = try!(fst::Map::from_path(fst));
    let by_id = try!(taxons.map(|taxons| taxon::read_taxa_file(taxons))
                           .map(|res| res.map(Some)).unwrap_or(Ok(None)))
        .map(|taxons| taxon::TaxonList::new(taxons));
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let line = try!(line);
        if line.starts_with('>') {
            println!("{}", line);
        } else if let Some(lca) = fst.get(&line) {
            if with_input {
                print!("{},", line);
            }
            if let Some(ref by_id) = by_id {

                let taxon = try!(by_id.get(lca as usize)
                                      .ok_or("LCA taxon id not in taxon list. Check compatibility with index."));
                println!("{},{},{}", taxon.id, taxon.name, taxon.rank);
            } else {
                println!("{}", lca);
            }
        }
    }
    Ok(())
}

fn prot2kmer2lca(fst: &str, k: &str) -> Result<()> {
    let map = try!(fst::Map::from_path(fst));
    let mut writer = fasta::Writer::new(io::stdout(), " ", false);
    let k = try!(k.parse::<usize>());

    for prot in fasta::Reader::new(io::stdin(), None, true).records() {
        let prot = try!(prot);

        if prot.sequence[0].len() < k {
            continue
        }

        let lcas = (0..(prot.sequence[0].len() - k + 1))
            .map(|i| &prot.sequence[0][i..i + k])
            .filter_map(|kmer| map.get(kmer))
            .map(|lca| lca.to_string())
            .collect::<Vec<_>>();

        if ! lcas.is_empty() {
            try!(writer.write_record(fasta::Record {
                header: prot.header,
                sequence: lcas
            }));
        }
    }

    Ok(())
}

fn taxa2agg(taxons: &str, method: &str, aggregation: &str, delimiter: Option<&str>, ranked_only: bool, factor: &str, scored: bool, lower_bound: &str) -> Result<()> {
    // Parsing the Taxa file
    let taxons = try!(taxon::read_taxa_file(taxons).map_err(|err| err.to_string()));

    // Parsing the factor
    let factor = try!(factor.parse::<f32>().map_err(|err| err.to_string()));

    // Parsing the factor
    let lower_bound = try!(lower_bound.parse::<f32>().map_err(|err| err.to_string()));

    // Parsing the delimiter regex
    let delimiter = try!(regex::Regex::new(delimiter.unwrap_or(r"\s+"))
                                      .map(Some)
                                      .map_err(|err| err.to_string()));

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    let aggregator: Result<Box<agg::Aggregator>> = match (method, aggregation) {
        ("RMQ",  "MRTL") => Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id))),
        ("RMQ",  "LCA*") => Ok(Box::new(rmq::lca::LCACalculator::new(tree))),
        ("RMQ",  "both") => {
            writeln!(&mut io::stderr(), "Warning: this is a hybrid between LCA/MRTL, not LCA*/MRTL").unwrap();
            Ok(Box::new(rmq::mix::MixCalculator::new(tree, factor)))
        },
        ("tree", "LCA*") => Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id))),
        ("tree", "both") => Ok(Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, factor))),
        _                => Err(ErrorKind::InvalidInvocation(format!("{} and {} cannot be combined", method, aggregation)).into())
    };
    let aggregator = try!(aggregator);

    fn with_score(pair: &String) -> Result<(TaxonId, f32)> {
        let split = pair.split('=').collect::<Vec<_>>();
        if split.len() != 2 { try!(Err("Taxon without score")); }
        Ok((try!(split[0].parse::<TaxonId>()), try!(split[1].parse::<f32>())))
    }

    fn not_scored(tid: &String) -> Result<(TaxonId, f32)> {
        Ok((try!(tid.parse::<TaxonId>()), 1.0))
    }

    let parser = if scored { with_score } else { not_scored };

    let mut writer = fasta::Writer::new(io::stdout(), ",", false);

    // Iterate over each read
    for record in fasta::Reader::new(io::stdin(), delimiter, true).records() {
        // Parse the sequence of LCA's
        let record = try!(record.map_err(|err| err.to_string()));
        let taxons = record.sequence.iter()
                                    .map(parser)
                                    .collect::<Result<Vec<(TaxonId, f32)>>>();
        let taxons = try!(taxons.map_err(|err| format!("Error reading taxons ({:?}): {}", record, err)));

        // Create a frequency table of taxons for this read (taking into account the lower bound)
        let counts = agg::count(taxons.into_iter());
        let counts = agg::filter(counts, lower_bound);

        // If we don't have a consensus taxon, leave out the read (i.e. consider undetected)
        if !counts.is_empty() {
            let aggregate = try!(aggregator.aggregate(&counts).map_err(|err| err.to_string()));
            let taxon = by_id.get(snapping[aggregate].unwrap()).unwrap();
            try!(writer.write_record(fasta::Record {
                header: record.header,
                sequence: vec![taxon.id.to_string(), taxon.name.to_string(), taxon.rank.to_string()],
            }).map_err(|err| err.to_string()));
        }
    }
    Ok(())
}

