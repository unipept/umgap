
use std::io;
use std::io::Write;
use std::io::BufRead;
use std::borrow::Cow;
use std::fs;
use std::collections::HashSet;
use std::collections::HashMap;
use std::ops;
use std::cmp;

#[macro_use(clap_app, crate_version, crate_authors)]
extern crate clap;

extern crate fst;

extern crate regex;

extern crate csv;

#[macro_use(quick_main)]
extern crate error_chain;

#[macro_use(json, json_internal)]
extern crate serde_json;
use serde_json::value;

extern crate either;
use either::Either;

extern crate umgap;
use umgap::dna::Strand;
use umgap::dna::translation::TranslationTable;
use umgap::io::fasta;
use umgap::io::fastq;
use umgap::taxon;
use umgap::errors::{Result, ErrorKind};
use umgap::taxon::TaxonId;
use umgap::agg;
use umgap::rmq;
use umgap::tree;
use umgap::utils;

quick_main!(|| -> Result<()> {

    // TODO https://docs.rs/clap/2.31.2/clap/macro.arg_enum.html
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
                "Adds a reading frame (1, 2, 3, 1R, 2R or 3R). (Default: only \
                 1)")
            (@arg append_name: -n --name
                "Append a bar (|) and the name of the frame to the fasta \
                header.")
            (@arg table: -t --table [INT]
                "Translation table to use (default: 1)")
            (@arg show_table: -s --("show-table")
                "Print the selected table and exit")
        )
        (@subcommand pept2lca =>
            (about: "Looks up each line of input in a given FST index and \
                     outputs the result. Lines starting with a '>' are \
                     copied. Lines for which no mapping is found are ignored.")
            (@arg one_on_one: -o --("one-on-one")
                "Map unknown sequences to 0 instead of ignoring them")
            (@arg fst_index: +required "An FST to query")
        )
        (@subcommand prot2kmer2lca =>
            (about: "Reads all the records in a specified FASTA file and \
                     queries the k-mers in an FST for the LCA's.")
            (@arg length: -k --length <INT>
                "The length of the k-mers in the FST")
            (@arg one_on_one: -o --("one-on-one")
                "Map unknown sequences to 0 instead of ignoring them")
            (@arg fst_index: +required "An FST to query")
        )
        (@subcommand taxa2agg =>
            (about: "Aggregates taxa to a single taxon.")
            (@arg scored: -s --scored
                "Each taxon is followed by a score between 0 and 1")
            (@arg ranked: -r --ranked
                "Restrict to taxa with a taxonomic rank")
            (@arg method: -m --method [STR]
                "The method to use for aggregation (RMQ or tree)")
            (@arg aggregate: -a --aggregate [STR]
                "The aggregation to use (LCA*, MRTL or hybrid)")
            (@arg factor: -f --factor [RATIO]
                "The factor for the hybrid aggregation, from 0.0 (MRTL) to \
                 1.0 (LCA*)")
            (@arg lower_bound: -l --("lower-bound") [INT]
                "The smallest input frequency for a taxon to be included in \
                 the aggregation")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
        (@subcommand prot2pept =>
            (about: "Splits each protein sequence in a FASTA format into a \
                     list of (tryptic) peptides.")
            (@arg pattern: -p --pattern [REGEX]
                "The cleavage-pattern (regex), i.e. the pattern after which \
                 the next peptide will be cleaved (default: ([KR])([^P]) for \
                 tryptic peptides).")
        )
        (@subcommand prot2kmer =>
            (about: "Splits each protein sequence in a FASTA format into a \
                     list of kmers.")
            (@arg length: -k --length <INT>
                "The K in K-mers")
        )
        (@subcommand uniq =>
            (about: "Concatenates the data strings of all consecutive FASTA \
                     entries with the same header.")
            (@arg output: -s --separator [INT]
                "Separator between output items (default the empty string)")
            (@arg input: -i --("input-separator") [INT]
                "Separator regex input items (default same as separator)")
            (@arg keep: -k --keep
                "Keep newline in the input sequence")
            (@arg wrap: -w --wrap
                "Wrap the output sequences")
        )
        (@subcommand filter =>
            (about: "Filter peptides in a FASTA format based on specific \
                     criteria.")
            (@arg min_length: -m --minlen [INT]
                "Only retain tryptic peptides that have at least min \
                 (default: 5) amino acids.")
            (@arg max_length: -M --maxlen [INT]
                "Only retain tryptic peptides that have at most max \
                 (default: 50) amino acids.")
            (@arg contains: -c --contains [STRING]
                "The letters that a sequence must contain.")
            (@arg lacks: -l --lacks [STRING]
                "The letters that a sequence cannot contain.")
        )
        (@subcommand fastq2fasta =>
            (about: "Interleaves a number of FASTQ files into a single FASTA \
                     output.")
            (@arg input: ... #{1,10} "The input files")
        )
        (@subcommand buildindex =>
            (about: "Write an FST index of stdin on stdout.")
        )
        (@subcommand snaprank =>
            (about: "Snap taxa to a specified rank.")
            (@arg rank: -r --rank [RANK]
                "The rank (default: species) to show.")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
        (@subcommand jsontree =>
            (about: "Aggregates taxa to a JSON tree for usage in the \
                     unipept visualisations.")
            (@arg ranked: -r --ranked
                "Restrict to taxa with a taxonomic rank")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
        (@subcommand seedextend =>
            (about: "Seed and extend.")
            (@arg min_seed_size: -s --("min-seed-size") [INT]
                "The minimum length of equal taxa to count as seed \
                 (default 4)")
            (@arg max_gap_size: -g --("max-gap-size") [INT]
                "The maximum length of a gap between seeds in an extension \
                 (default 0)")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
        (@subcommand report =>
            (about: "Count and report on a list of taxon ids")
            (@arg rank: -r --rank [RANK]
                "The rank (default: species) to show.")
            (@arg taxon_file: +required "The NCBI taxonomy tsv-file")
        )
        (@subcommand bestof =>
            (about: "Pick the frame with the most none-root hits")
            (@arg frames: -f --frames [INT]
                "The number of frames of which to pick the best")
        )
    ).get_matches();

    match matches.subcommand() {
        ("translate", Some(matches)) => translate(
            matches.is_present("methionine"),
            matches.value_of("table").unwrap_or("1"),
            matches.is_present("show_table"),
            matches.is_present("append_name"),
            matches.values_of("frame").map(Iterator::collect).unwrap_or(
                if matches.is_present("all_frames") { vec!["1", "2", "3", "1R", "2R", "3R"] } else { vec!["1"] }
            )),
        ("pept2lca", Some(matches)) => pept2lca(
            matches.value_of("fst_index").unwrap(), // required so safe
            matches.is_present("one_on_one")),
        ("prot2kmer2lca", Some(matches)) => prot2kmer2lca(
            matches.value_of("fst_index").unwrap(), // required so safe
            matches.value_of("length").unwrap(),
            matches.is_present("one_on_one")),
        ("taxa2agg", Some(matches)) => taxa2agg(
            matches.value_of("taxon_file").unwrap(), // required so safe
            matches.value_of("method").unwrap_or("RMQ"),
            matches.value_of("aggregate").unwrap_or("LCA*"),
            matches.is_present("ranked"),
            matches.value_of("factor").unwrap_or("0"),
            matches.is_present("scored"),
            matches.value_of("lower_bound").unwrap_or("0")),
        ("prot2pept", Some(matches)) => prot2pept(
            matches.value_of("pattern").unwrap_or("([KR])([^P])")),
        ("prot2kmer", Some(matches)) => prot2kmer(
            matches.value_of("length").unwrap()),
        ("uniq", Some(matches)) => uniq(
            matches.value_of("output").unwrap_or(""),
            matches.value_of("input"),
            matches.is_present("keep"),
            matches.is_present("wrap")),
        ("filter", Some(matches)) => filter(
            matches.value_of("min_length").unwrap_or("5"),
            matches.value_of("max_length").unwrap_or("50"),
            matches.value_of("contains").unwrap_or(""),
            matches.value_of("lacks").unwrap_or("")),
        ("fastq2fasta", Some(matches)) => fastq2fasta(
            matches.values_of("input").unwrap().collect()), // required so safe
        ("buildindex", Some(_)) => buildindex(),
        ("snaprank", Some(matches)) => snaprank(
            matches.value_of("taxon_file").unwrap(), // required so safe
            matches.value_of("rank").unwrap_or("species")),
        ("jsontree", Some(matches)) => jsontree(
            matches.value_of("taxon_file").unwrap(), // required so safe
            matches.is_present("ranked")),
        ("seedextend", Some(matches)) => seedextend(
            matches.value_of("min_seed_size").unwrap_or("4"),
            matches.value_of("max_gap_size").unwrap_or("0"),
            matches.value_of("taxon_file").unwrap()), // required so safe
        ("report", Some(matches)) => report(
            matches.value_of("taxon_file").unwrap(), // required so safe
            matches.value_of("rank").unwrap_or("species")),
        ("bestof", Some(matches)) => bestof(
            matches.value_of("frames").unwrap_or("6")),
        _  => { println!("{}", matches.usage()); Ok(()) }
    }
});

fn translate(methionine: bool, table: &str, show_table: bool, append_name: bool, frames: Vec<&str>) -> Result<()> {
    // Parsing the table
    let table = table.parse::<&TranslationTable>()?;

    // Split on show_tables
    if show_table {
        table.print();
    } else {
        let mut writer = fasta::Writer::new(io::stdout(), "", false);

        // Parsing the frames
        let frames = frames.iter().map(|&frame| match frame {
            "1"  => Ok((frame, 1, false)),
            "2"  => Ok((frame, 2, false)),
            "3"  => Ok((frame, 3, false)),
            "1R" => Ok((frame, 1, true)),
            "2R" => Ok((frame, 2, true)),
            "3R" => Ok((frame, 3, true)),
            _    => Err(ErrorKind::InvalidInvocation(format!("{} is not a frame", frame)).into())
        }).collect::<Result<Vec<(&str, usize, bool)>>>()?;

        for record in fasta::Reader::new(io::stdin(), None, true).records() {
            let fasta::Record { header, sequence } = record?;

            let forward = Strand::from(sequence[0].as_bytes());
            let reverse = forward.reversed();
            for &(name, frame, reversed) in &frames {
                let strand = if reversed { &reverse } else { &forward };
                writer.write_record(fasta::Record {
                    header: if !append_name { header.clone() } else { header.clone() + "|" + name },
                    sequence: vec![String::from_utf8(table.translate_frame(methionine, strand.frame(frame))).unwrap()]
                })?;
            }
        }
    }
    Ok(())
}

fn pept2lca(fst: &str, one_on_one: bool) -> Result<()> {
    let fst = fst::Map::from_path(fst)?;
    let default = if one_on_one { Some(0) } else { None };
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(stdout, "{}", line)?;
        } else if let Some(lca) = fst.get(&line).map(Some).unwrap_or(default) {
            writeln!(stdout, "{}", lca)?;
        }
    }
    Ok(())
}

fn prot2kmer2lca(fst: &str, k: &str, one_on_one: bool) -> Result<()> {
    let map = fst::Map::from_path(fst)?;
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    let k = k.parse::<usize>()?;

    for prot in fasta::Reader::new(io::stdin(), None, true).records() {
        let prot = prot?;

        if prot.sequence[0].len() < k {
            continue
        }

        let lcas = {
            let kmers = (0..(prot.sequence[0].len() - k + 1)).map(|i| &prot.sequence[0][i..i + k]);
            if one_on_one {
                Either::Left(kmers.map(|kmer| map.get(kmer).unwrap_or(0)))
            } else {
                Either::Right(kmers.filter_map(|kmer| map.get(kmer)))
            }
        }.map(|lca| lca.to_string())
         .collect::<Vec<_>>();

        if !lcas.is_empty() {
            writer.write_record(fasta::Record {
                header: prot.header,
                sequence: lcas
            })?;
        }
    }

    Ok(())
}

fn taxa2agg(taxons: &str, method: &str, aggregation: &str, ranked_only: bool, factor: &str, scored: bool, lower_bound: &str) -> Result<()> {
    // Parsing the Taxa file
    let taxons = taxon::read_taxa_file(taxons)?;

    // Parsing the factor
    let factor = factor.parse::<f32>()?;

    // Parsing the factor
    let lower_bound = lower_bound.parse::<f32>()?;

    // Parsing the delimiter regex
    let delimiter = Some(regex::Regex::new("\n").unwrap());

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    let aggregator: Result<Box<agg::Aggregator>> = match (method, aggregation) {
        ("RMQ",  "MRTL") => Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id))),
        ("RMQ",  "LCA*") => Ok(Box::new(rmq::lca::LCACalculator::new(tree))),
        ("RMQ",  "hybrid") => {
            writeln!(&mut io::stderr(), "Warning: this is a hybrid between LCA/MRTL, not LCA*/MRTL").unwrap();
            Ok(Box::new(rmq::mix::MixCalculator::new(tree, factor)))
        },
        ("tree", "LCA*") => Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id))),
        ("tree", "hybrid") => Ok(Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, factor))),
        _                => Err(ErrorKind::InvalidInvocation(format!("{} and {} cannot be combined", method, aggregation)).into())
    };
    let aggregator = aggregator?;

    fn with_score(pair: &String) -> Result<(TaxonId, f32)> {
        let split = pair.split('=').collect::<Vec<_>>();
        if split.len() != 2 { Err("Taxon without score")?; }
        Ok((split[0].parse::<TaxonId>()?, split[1].parse::<f32>()?))
    }

    fn not_scored(tid: &String) -> Result<(TaxonId, f32)> {
        Ok((tid.parse::<TaxonId>()?, 1.0))
    }

    let parser = if scored { with_score } else { not_scored };

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    // Iterate over each read
    for record in fasta::Reader::new(io::stdin(), delimiter, true).records() {
        // Parse the sequence of LCA's
        let record = record?;
        let taxons = record.sequence.iter()
                           .map(parser)
                           .collect::<Result<Vec<(TaxonId, f32)>>>()?;

        // Create a frequency table of taxons for this read (taking into account the lower bound)
        let counts = agg::count(taxons.into_iter());
        let counts = agg::filter(counts, lower_bound);

        writer.write_record(fasta::Record {
            header: record.header,
            sequence: if counts.is_empty() { vec![] } else {
                let aggregate = aggregator.aggregate(&counts)?;
                vec![snapping[aggregate].unwrap().to_string()]
            }
        })?;
    }
    Ok(())
}

fn prot2pept(pattern: &str) -> Result<()> {
    let pattern = regex::Regex::new(pattern)?;

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = record?;

        // We will run the regex replacement twice, since a letter can be
        // matched twice (i.e. once after and once before the split).
        let first_run = pattern.replace_all(&sequence[0], "$1\n$2");

        writer.write_record(fasta::Record {
            header: header,
            sequence: pattern.replace_all(&first_run, "$1\n$2")
                             .replace("*", "\n")
                             .lines()
                             .filter(|x| !x.is_empty())
                             .map(ToOwned::to_owned)
                             .collect(),
        })?;
    }
    Ok(())
}

fn prot2kmer(k: &str) -> Result<()> {
    let k = k.parse::<usize>()?;

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = record?;
        if sequence[0].len() < k { continue }
        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence[0].as_bytes().windows(k)
                                 .map(String::from_utf8_lossy).map(Cow::into_owned)
                                 .collect(),
        })?;
    }
    Ok(())
}

fn filter(min_length: &str, max_length: &str, contains: &str, lacks: &str) -> Result<()> {
    let min      = min_length.parse::<usize>()?;
    let max      = max_length.parse::<usize>()?;
    let contains = contains.chars().collect::<HashSet<char>>();
    let lacks    = lacks.chars().collect::<HashSet<char>>();

    // Each peptide/nucleotide sequence is assumed to be on its own line
    let delimiter = regex::Regex::new(r"\n").map(Some)?;

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), delimiter, false).records() {
        let fasta::Record { header, sequence } = record?;

        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence.into_iter()
                              .filter(|seq| {
                                  let length = seq.len();
                                  length >= min && length <= max
                              })
                              .filter(|seq| {
                                  let set = seq.chars().collect::<HashSet<char>>();
                                  contains.intersection(&set).count() == contains.len()
                                    && lacks.intersection(&set).count() == 0
                              })
                              .collect(),
        })?;
    }
    Ok(())
}

fn uniq(separator: &str, input_separator: Option<&str>, keep: bool, wrap: bool) -> Result<()> {
    // Parsing the input separator regex
    let input_separator = regex::Regex::new(input_separator.unwrap_or(&regex::escape(separator)))?;

    let mut last   = None::<fasta::Record>;
    let mut writer = fasta::Writer::new(io::stdout(), separator, wrap);
    for record in fasta::Reader::new(io::stdin(), Some(input_separator), !keep).records() {
        let record = record?;
        if let Some(ref mut rec) = last {
            if rec.header == record.header {
                rec.sequence.extend(record.sequence);
            } else {
                writer.write_record_ref(rec)?;
                *rec = record;
            }
        } else {
            last = Some(record);
        }
    }
    if let Some(rec) = last {
        writer.write_record(rec)?;
    }
    Ok(())
}

fn fastq2fasta(input: Vec<&str>) -> Result<()> {
    let handles = input.iter()
                       .map(fs::File::open)
                       .collect::<io::Result<Vec<fs::File>>>()?;
    let readers = handles.iter()
                         .map(fastq::Reader::new)
                         .map(fastq::Reader::records)
                         .collect();
    let mut writer = fasta::Writer::new(io::stdout(), "", false);
    for recordzip in utils::Zip::new(readers) {
        for record in recordzip {
            let record = record?;
            writer.write_record(fasta::Record {
                header: record.header,
                sequence: vec![record.sequence],
            })?;
        }
    }
    Ok(())
}

fn buildindex() -> Result<()> {
    let mut reader = csv::Reader::from_reader(io::stdin())
                                 .has_headers(false)
                                 .delimiter(b'\t');

    let mut index = fst::MapBuilder::new(io::stdout())?;

    for record in reader.decode() {
        let (kmer, lca): (String, u64) = record?;
        index.insert(kmer, lca)?;
    }

    index.finish()?;

    Ok(())
}

fn snaprank(taxons: &str, rank: &str) -> Result<()> {
    let taxons = taxon::read_taxa_file(taxons)?;
    let rank = rank.parse::<taxon::Rank>()?;
    if rank == taxon::Rank::NoRank {
        return Err(ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.filter_ancestors(|tid|
        by_id.get(tid).map(|t| t.rank == rank).unwrap_or(false)
    );

    // Read and count taxon ranks
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(handle, "{}", line)?;
        } else {
            let taxon = line.parse::<taxon::TaxonId>()?;
            let snapped = snapping[taxon].unwrap_or(0);
            writeln!(handle, "{}", snapped)?;
        }
    }

    Ok(())
}

fn jsontree(taxons: &str, ranked_only: bool) -> Result<()> {
    let taxons = taxon::read_taxa_file(taxons)?;

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, ranked_only);

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = line?.parse::<taxon::TaxonId>()?;
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    // Recursive json transformation
    fn to_json(node: &tree::tree::Tree<usize>, aggnode: &tree::tree::Tree<usize>, by_id: &taxon::TaxonList) -> value::Value {
        let root = by_id.get(node.root).unwrap();
        json!({
            "name": root.name,
            "id": node.root,
            "data": {
                "count": aggnode.value,
                "valid_taxon": if root.valid { "1" } else { "0" },
                "rank": format!("{}", root.rank),
                "self_count": node.value
            },
            "children": node.children.iter().zip(aggnode.children.iter())
                                     .map(|(n, s)| to_json(n, s, by_id))
                                     .collect::<Vec<_>>()
        })
    }

    let tree = tree::tree::Tree::new(1, &by_id.ancestry(), &counts)?;
    let aggtree = tree.aggregate(&ops::Add::add);
    print!("{}", to_json(&tree, &aggtree, &by_id));

    Ok(())
}

fn seedextend(min_seed_size: &str, max_gap_size: &str, _taxon_file: &str) -> Result<()> {
    let min_seed_size = min_seed_size.parse::<usize>()?;
    let max_gap_size = max_gap_size.parse::<usize>()?;
    let _gap_penalty = 0.5; // TODO: move down
    // let _taxa = taxon::read_taxa_file(taxon_file)?; // TODO use it

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    let separator = Some(regex::Regex::new("\n").unwrap());
    for record in fasta::Reader::new(io::stdin(), separator, false).records() {
        let record = record?;
        let taxons = record.sequence.iter()
                                    .map(|s| s.parse::<TaxonId>())
                                    .collect::<std::result::Result<Vec<TaxonId>,_>>()?;
        if taxons.len() == 0 {
            writer.write_record(fasta::Record {
                header: record.header,
                sequence: vec![]
            })?;
            continue;
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

fn report(taxons: &str, rank: &str) -> Result<()> {
    let taxons = taxon::read_taxa_file(taxons)?;
    let rank = rank.parse::<taxon::Rank>()?;
    if rank == taxon::Rank::NoRank {
        return Err(ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }

    // Parsing the taxons
    let tree     = taxon::TaxonTree::new(&taxons);
    let by_id    = taxon::TaxonList::new(taxons);
    let snapping = tree.filter_ancestors(|tid|
        by_id.get(tid).map(|t| t.rank == rank).unwrap_or(false)
    );

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = line?.parse::<taxon::TaxonId>()?;
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    let mut counts = counts
        .into_iter()
        .map(|(taxon, count)| (count, taxon))
        .collect::<Vec<(taxon::TaxonId, usize)>>();
    counts.sort();

    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    for (count, taxon) in counts {
        let taxon = by_id.get(taxon)
                         .ok_or("LCA taxon id not in taxon list. Check compatibility with index.")?;
        writeln!(stdout, "{},{},{}", count, taxon.id, taxon.name)?;
    }

    Ok(())
}

fn bestof(frames: &str) -> Result<()> {
    // Parsing the number of frames
    let frames = frames.parse::<usize>()?;

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    let delimiter = Some(regex::Regex::new(r"\n").unwrap());

    // Combine frames and process them
    let mut chunk = Vec::with_capacity(frames);
    for record in fasta::Reader::new(io::stdin(), delimiter, false).records() {
        let record = record?;
        if chunk.len() < frames - 1 {
            chunk.push(record);
        } else {
            // process chunk
            writer.write_record_ref(chunk.iter().max_by_key(|&rec| {
                rec.sequence.iter()
                   .map(|tid| tid.parse::<TaxonId>().unwrap_or(0))
                   .filter(|&s| s != 0 && s != 1)
                   .count()
            }).unwrap())?;
            chunk.clear();
        }
    }
    Ok(())
}
