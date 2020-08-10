use std::borrow::Cow;
use std::cmp;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::ops;

use fst;
use fst::Streamer;

use csv;

#[macro_use(quick_main)]
extern crate error_chain;

#[macro_use(json)]
extern crate serde_json;
use serde_json::value;

use structopt::StructOpt;

use umgap::agg;
use umgap::agg::Aggregator;
use umgap::args;
use umgap::commands;
use umgap::errors::{ErrorKind, Result};
use umgap::io::fasta;
use umgap::io::fastq;
use umgap::rank;
use umgap::taxon;
use umgap::taxon::TaxonId;
use umgap::tree;
use umgap::utils;

quick_main!(|| -> Result<()> {
    match args::Opt::from_args() {
        args::Opt::Translate(args) => commands::translate::translate(args),
        args::Opt::PeptToLca(args) => commands::pept2lca::pept2lca(args),
        args::Opt::ProtToKmerToLca(args) => commands::prot2kmer2lca::prot2kmer2lca(args),
        args::Opt::ProtToTrypToLca(args) => commands::prot2tryp2lca::prot2tryp2lca(args),
        args::Opt::TaxaToAgg(args) => commands::taxa2agg::taxa2agg(args),
        args::Opt::ProtToPept(args) => commands::prot2pept::prot2pept(args),
        args::Opt::ProtToKmer(args) => prot2kmer(args),
        args::Opt::Filter(args) => filter(args),
        args::Opt::Uniq(args) => uniq(args),
        args::Opt::FastqToFasta(args) => fastq2fasta(args),
        args::Opt::Taxonomy(args) => taxonomy(args),
        args::Opt::SnapTaxon(args) => snaptaxon(args),
        args::Opt::JsonTree(args) => jsontree(args),
        args::Opt::SeedExtend(args) => seedextend(args),
        args::Opt::Report(args) => commands::report::report(args),
        args::Opt::BestOf(args) => commands::bestof::bestof(args),
        args::Opt::PrintIndex(args) => printindex(args),
        args::Opt::SplitKmers(args) => splitkmers(args),
        args::Opt::JoinKmers(args) => joinkmers(args),
        args::Opt::BuildIndex => buildindex(),
        args::Opt::CountRecords => countrecords(),
        args::Opt::Visualize(args) => visualize(args),
    }
});

fn prot2kmer(args: args::ProtToKmer) -> Result<()> {
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), true).records() {
        let fasta::Record { header, sequence } = record?;
        if sequence[0].len() < args.length {
            continue;
        }
        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence[0]
                .as_bytes()
                .windows(args.length)
                .map(String::from_utf8_lossy)
                .map(Cow::into_owned)
                .collect(),
        })?;
    }
    Ok(())
}

fn filter(args: args::Filter) -> Result<()> {
    let contains = args.contains.chars().collect::<HashSet<char>>();
    let lacks = args.lacks.chars().collect::<HashSet<char>>();

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let fasta::Record { header, sequence } = record?;

        writer.write_record(fasta::Record {
            header: header,
            sequence: sequence
                .into_iter()
                .filter(|seq| {
                    let length = seq.len();
                    length >= args.min_length && length <= args.max_length
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

fn uniq(args: args::Uniq) -> Result<()> {
    let mut last = None::<fasta::Record>;
    let mut writer = fasta::Writer::new(io::stdout(), &args.separator, args.wrap);
    for record in fasta::Reader::new(io::stdin(), false).records() {
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

fn fastq2fasta(args: args::FastqToFasta) -> Result<()> {
    let handles = args
        .input
        .iter()
        .map(fs::File::open)
        .collect::<io::Result<Vec<fs::File>>>()?;
    let readers = handles
        .iter()
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

fn taxonomy(args: args::Taxonomy) -> Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    let by_id = taxon::TaxonList::new(taxons);

    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if !args.no_header {
        write!(handle, "taxon_id,taxon_name,taxon_rank")?;
        if args.all_ranks {
            for rank in rank::Rank::ranks() {
                let rank_name = rank.to_string().replace(" ", "_");
                write!(handle, ",{}_id,{}_name", rank_name, rank_name)?;
            }
        }
        write!(handle, "\n")?;
    }

    for line in stdin.lock().lines() {
        let line = line?;
        if line.starts_with('>') {
            writeln!(handle, "{}", line)?;
        } else {
            let id = line.parse::<taxon::TaxonId>()?;
            // Map to root if id not found
            let taxon = by_id.get_or_unknown(id)?;
            write!(handle, "{},{},{}", taxon.id, taxon.name, taxon.rank)?;
            if args.all_ranks {
                let lineage = by_id.lineage(id)?;
                for rank in rank::Rank::ranks() {
                    if let Some(l_taxon) = &lineage[rank] {
                        write!(handle, ",{},{}", l_taxon.id, l_taxon.name)?;
                    } else {
                        write!(handle, ",,")?;
                    }
                }
            }
            write!(handle, "\n")?;
        }
    }
    Ok(())
}

fn snaptaxon(args: args::SnapTaxon) -> Result<()> {
    let taxons = taxon::read_taxa_file(&args.taxon_file)?;
    if args.rank.map(|r| r == rank::Rank::NoRank).unwrap_or(false) {
        return Err(ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
    }

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = tree.filter_ancestors(|tid| {
        args.taxons.contains(&tid)
            || by_id
                .get(tid)
                .map(|t| {
                    (args.invalid || t.valid) && args.rank.map(|r| t.rank == r).unwrap_or(false)
                })
                .unwrap_or(false)
    });

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

fn jsontree(args: args::JsonTree) -> Result<()> {
    let taxons = taxon::read_taxa_file(args.taxon_file)?;

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, args.ranked_only);

    // Read and count taxon ranks
    let mut counts = HashMap::new();
    let stdin = io::stdin();
    for line in stdin.lock().lines() {
        let taxon = line?.parse::<taxon::TaxonId>()?;
        *counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
    }

    // Recursive json transformation
    fn to_json(
        node: &tree::tree::Tree<usize>,
        aggnode: &tree::tree::Tree<usize>,
        by_id: &taxon::TaxonList,
        freq: usize,
    ) -> value::Value {
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
                                     .filter(|&(_n, s)| s.value > freq)
                                     .map(|(n, s)| to_json(n, s, by_id, freq))
                                     .collect::<Vec<_>>()
        })
    }

    let tree = tree::tree::Tree::new(1, &by_id.ancestry(), &counts)?;
    let aggtree = tree.aggregate(&ops::Add::add);
    print!("{}", to_json(&tree, &aggtree, &by_id, args.min_frequency));

    Ok(())
}

fn seedextend(args: args::SeedExtend) -> Result<()> {
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    let by_id = if let Some(ref tf) = args.ranked {
        let taxa = taxon::read_taxa_file(tf)?;
        Some(taxon::TaxonList::new_with_unknown(taxa, true))
    } else {
        None
    };

    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        let mut taxons = record
            .sequence
            .iter()
            .map(|s| s.parse::<TaxonId>())
            .collect::<std::result::Result<Vec<TaxonId>, _>>()?;
        taxons.push(0);

        let mut seeds = Vec::new();
        let mut start = 0;
        let mut end = 1;
        let mut last_tid = taxons[start];
        let mut same_tid = 1;
        let mut same_max = 1;
        while end < taxons.len() {
            // same tid as last, add to seed.
            if last_tid == taxons[end] {
                same_tid += 1;
                end += 1;
                continue;
            }

            // our gap just became to big
            if last_tid == 0 && same_tid > args.max_gap_size {
                // add extended seed
                if same_max >= args.min_seed_size {
                    seeds.push((start, end - same_tid))
                }
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
            if last_tid != 0 {
                same_max = cmp::max(same_max, same_tid);
            }
            last_tid = taxons[end];
            same_tid = 1;
            end += 1;
        }
        if same_max >= args.min_seed_size {
            if last_tid == 0 {
                end -= same_tid
            }
            seeds.push((start, end))
        }

        if let Some(ref by_id) = by_id {
            seeds = seeds
                .into_iter()
                .max_by_key(|(s, e)| {
                    taxons
                        .iter()
                        .skip(*s)
                        .take(e - s)
                        .map(|t| by_id.score(*t).unwrap_or(args.penalty))
                        .sum::<usize>()
                })
                .into_iter()
                .collect::<Vec<(usize, usize)>>();
        }

        // write it
        writer
            .write_record(fasta::Record {
                header: record.header,
                sequence: seeds
                    .into_iter()
                    .flat_map(|(start, end)| taxons[start..end].into_iter())
                    .map(|t| t.to_string())
                    .collect(),
            })
            .map_err(|err| err.to_string())?;
    }
    Ok(())
}

fn printindex(args: args::PrintIndex) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    let index = unsafe { fst::Map::from_path(args.fst_file) }?;
    let mut stream = index.stream();

    while let Some((k, v)) = stream.next() {
        writer.serialize((String::from_utf8_lossy(k), v))?;
    }

    Ok(())
}

fn splitkmers(args: args::SplitKmers) -> Result<()> {
    let byte = args.prefix.as_bytes().first();

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());

    for record in reader.deserialize() {
        let (tid, sequence): (TaxonId, String) = record?;
        if sequence.len() < args.length {
            continue;
        }
        for kmer in sequence.as_bytes().windows(args.length) {
            if byte.is_none() {
                writer.serialize((String::from_utf8_lossy(kmer), tid))?;
            } else if *byte.unwrap() == kmer[0] {
                writer.serialize((String::from_utf8_lossy(&kmer[1..]), tid))?;
            }
        }
    }

    Ok(())
}

fn buildindex() -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());

    let mut index = fst::MapBuilder::new(io::stdout())?;

    for record in reader.deserialize() {
        let (kmer, lca): (String, u64) = record?;
        index.insert(kmer, lca)?;
    }

    index.finish()?;

    Ok(())
}

fn joinkmers(args: args::JoinKmers) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let taxons = taxon::read_taxa_file(args.taxon_file)?;

    // Parsing the Taxa file
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let ranksnapping = tree.snapping(&by_id, true);
    let validsnapping = tree.snapping(&by_id, false);
    let aggregator = tree::mix::MixCalculator::new(tree.root, &by_id, 0.95);

    let mut emit = |kmer: &str, tids: Vec<(TaxonId, f32)>| {
        let counts = agg::count(tids.into_iter());
        if let Ok(aggregate) = aggregator.aggregate(&counts) {
            let taxon = ranksnapping[aggregate].unwrap();
            let rank = by_id.get_or_unknown(taxon).unwrap().rank;
            write!(handle, "{}\t{}\t{}\n", kmer, taxon, rank)
        } else {
            Ok(())
        }
    };

    // Iterate over records and emit groups
    let mut current_kmer: Option<String> = Option::None;
    let mut current_tids = vec![];
    for record in reader.deserialize() {
        let (kmer, tid): (String, TaxonId) = record?;
        if let Some(c) = current_kmer {
            if c != kmer {
                emit(&c, current_tids)?;
                current_tids = vec![];
            }
        } else {
            current_tids = vec![];
        }
        current_kmer = Some(kmer);
        if let Some(validancestor) = validsnapping[tid] {
            current_tids.push((validancestor, 1.0));
        }
    }
    if let Some(c) = current_kmer {
        emit(&c, current_tids)?;
    }

    Ok(())
}

fn countrecords() -> Result<()> {
    let mut records = 0;
    let mut sequences = 0;
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        records += 1;
        for seq in record.sequence {
            if seq.len() > 0 {
                sequences += 1;
            }
        }
    }
    println!("Records: {}", records);
    println!("Sequence items: {}", sequences);
    Ok(())
}

fn visualize(args: args::Visualize) -> Result<()> {
    let mut taxa = HashMap::new();
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        let taxon = record.sequence[0].parse::<TaxonId>()?;
        *taxa.entry(taxon).or_insert(0) += 1;
    }

    let json = json!({
        "counts": taxa,
        "link": args.url.to_string(),
    });

    let client = reqwest::blocking::Client::new();
    let res = client
        .post("http://api.unipept.ugent.be/api/v1/taxa2tree")
        .json(&json)
        .send()
        .map_err(|err| err.to_string())?;

    print!("{}", res.text().map_err(|err| err.to_string())?);
    Ok(())
}
