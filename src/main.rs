#![cfg_attr(rustfmt, rustfmt_skip)]

use std::io;
use std::io::Write;
use std::io::BufRead;
use std::borrow::Cow;
use std::fs;
use std::collections::HashSet;
use std::collections::HashMap;
use std::ops;
use std::cmp;

extern crate clap;

extern crate fst;
use fst::Streamer;

extern crate regex;

extern crate csv;

#[macro_use(quick_main)]
extern crate error_chain;

#[macro_use(json, json_internal)]
extern crate serde_json;
use serde_json::value;

extern crate either;
use either::Either;

extern crate structopt;
use structopt::StructOpt;

extern crate strum;

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
use umgap::args;
use umgap::rank;


quick_main!(|| -> Result<()> {
	match args::Opt::from_args() {
		args::Opt::Translate(args)       => translate(args),
		args::Opt::PeptToLca(args)       => pept2lca(args),
		args::Opt::ProtToKmerToLca(args) => prot2kmer2lca(args),
		args::Opt::TaxaToAgg(args)       => taxa2agg(args),
		args::Opt::ProtToPept(args)      => prot2pept(args),
		args::Opt::ProtToKmer(args)      => prot2kmer(args),
		args::Opt::Filter(args)          => filter(args),
		args::Opt::Uniq(args)            => uniq(args),
		args::Opt::FastqToFasta(args)    => fastq2fasta(args),
		args::Opt::Taxonomy(args)        => taxonomy(args),
		args::Opt::SnapRank(args)        => snaprank(args),
		args::Opt::JsonTree(args)        => jsontree(args),
		args::Opt::SeedExtend(args)      => seedextend(args),
		args::Opt::Report(args)          => report(args),
		args::Opt::BestOf(args)          => bestof(args),
		args::Opt::PrintIndex(args)      => printindex(args),
		args::Opt::BuildIndex            => buildindex(),
		args::Opt::CountRecords          => countrecords(),
	}
});

fn translate(args: args::Translate) -> Result<()> {
	// Parsing the table
	let table = args.table.parse::<&TranslationTable>()?;

	// Which frames TODO ugly
	let frames = if args.all_frames { vec![args::Frame::Forward1, args::Frame::Forward2, args::Frame::Forward3,
	                                       args::Frame::Reverse1, args::Frame::Reverse2, args::Frame::Reverse3] }
	                           else { args.frames };

	// Split on show_tables
	if args.show_table {
		table.print();
	} else {
		let mut writer = fasta::Writer::new(io::stdout(), "", false);

		// Parsing the frames
		let frames = frames.iter().map(|&frame| match frame {
		    args::Frame::Forward1 => (frame, 1, false),
		    args::Frame::Forward2 => (frame, 2, false),
		    args::Frame::Forward3 => (frame, 3, false),
		    args::Frame::Reverse1 => (frame, 1, true),
		    args::Frame::Reverse2 => (frame, 2, true),
		    args::Frame::Reverse3 => (frame, 3, true),
		}).collect::<Vec<(args::Frame, usize, bool)>>();

		for record in fasta::Reader::new(io::stdin(), true).records() {
			let fasta::Record { header, sequence } = record?;

			let forward = Strand::from(sequence[0].as_bytes());
			let reverse = forward.reversed();
			for &(name, frame, reversed) in &frames {
				let strand = if reversed { &reverse } else { &forward };
				writer.write_record(fasta::Record {
				    header: if !args.append_name { header.clone() } else { header.clone() + "|" + &name.to_string() },
				    sequence: vec![String::from_utf8(table.translate_frame(args.methionine, strand.frame(frame))).unwrap()]
				})?;
			}
		}
	}
	Ok(())
}

fn pept2lca(args: args::PeptToLca) -> Result<()> {
	let fst = unsafe { fst::Map::from_path(args.fst_file) }?;
	let default = if args.one_on_one { Some(0) } else { None };
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

fn prot2kmer2lca(args: args::ProtToKmerToLca) -> Result<()> {
	let map = unsafe { fst::Map::from_path(&args.fst_file) }?;
	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

	for prot in fasta::Reader::new(io::stdin(), true).records() {
		let prot = prot?;

		if prot.sequence[0].len() < args.length {
			continue;
		}

		let lcas = {
			let kmers = (0..(prot.sequence[0].len() - args.length + 1)).map(|i| &prot.sequence[0][i..i + args.length]);
			if args.one_on_one {
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

fn taxa2agg(args: args::TaxaToAgg) -> Result<()> {
	// Parsing the Taxa file
	let taxons = taxon::read_taxa_file(args.taxon_file)?;

	// Parsing the taxons
	let tree	 = taxon::TaxonTree::new(&taxons);
	let by_id	= taxon::TaxonList::new(taxons);
	let snapping = tree.snapping(&by_id, args.ranked_only);

	let aggregator: Result<Box<agg::Aggregator>> = match (args.method, args.strategy) {
		(args::Method::RangeMinimumQuery, args::Strategy::MaximumRootToLeafPath) => Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id))),
		(args::Method::RangeMinimumQuery, args::Strategy::LowestCommonAncestor) => Ok(Box::new(rmq::lca::LCACalculator::new(tree))),
		(args::Method::RangeMinimumQuery, args::Strategy::Hybrid) => {
			writeln!(&mut io::stderr(), "Warning: this is a hybrid between LCA/MRTL, not LCA*/MRTL").unwrap();
			Ok(Box::new(rmq::mix::MixCalculator::new(tree, args.factor)))
		},
		(args::Method::Tree, args::Strategy::LowestCommonAncestor) => Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id))),
		(args::Method::Tree, args::Strategy::Hybrid) => Ok(Box::new(tree::mix::MixCalculator::new(tree.root, &by_id, args.factor))),
		(m, s) => Err(ErrorKind::InvalidInvocation(format!("{:?} and {:?} cannot be combined", m, s)).into())
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

	let parser = if args.scored { with_score } else { not_scored };

	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

	// Iterate over each read
	for record in fasta::Reader::new(io::stdin(), false).records() {
		// Parse the sequence of LCA's
		let record = record?;
		let taxons = record.sequence.iter()
		                   .map(parser)
		                   .collect::<Result<Vec<(TaxonId, f32)>>>()?;

		// Create a frequency table of taxons for this read (taking into account the lower bound)
		let counts = agg::count(taxons.into_iter().filter(|&(tid, _)| tid != 0));
		let counts = agg::filter(counts, args.lower_bound);

		writer.write_record(fasta::Record {
			header: record.header,
			sequence: if counts.is_empty() { vec!["1".into()] } else {
				let aggregate = aggregator.aggregate(&counts)?;
				vec![snapping[aggregate].unwrap().to_string()]
			}
		})?;
	}
	Ok(())
}

fn prot2pept(args: args::ProtToPept) -> Result<()> {
	let pattern = regex::Regex::new(&args.pattern)?;

	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
	for record in fasta::Reader::new(io::stdin(), true).records() {
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

fn prot2kmer(args: args::ProtToKmer) -> Result<()> {
	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
	for record in fasta::Reader::new(io::stdin(), true).records() {
		let fasta::Record { header, sequence } = record?;
		if sequence[0].len() < args.length { continue }
		writer.write_record(fasta::Record {
			header: header,
			sequence: sequence[0].as_bytes().windows(args.length)
			                     .map(String::from_utf8_lossy).map(Cow::into_owned)
			                     .collect(),
		})?;
	}
	Ok(())
}

fn filter(args: args::Filter) -> Result<()> {
	let contains = args.contains.chars().collect::<HashSet<char>>();
	let lacks	= args.lacks.chars().collect::<HashSet<char>>();

	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);
	for record in fasta::Reader::new(io::stdin(), false).records() {
		let fasta::Record { header, sequence } = record?;

		writer.write_record(fasta::Record {
			header: header,
			sequence: sequence.into_iter()
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
	let mut last   = None::<fasta::Record>;
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
	let handles = args.input.iter()
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

fn taxonomy(args: args::Taxonomy) -> Result<()> {
	 let taxons = taxon::read_taxa_file(&args.taxon_file)?;
	 let by_id	= taxon::TaxonList::new(taxons);

	 let stdin = io::stdin();
	 let stdout = io::stdout();
	 let mut handle = stdout.lock();

	 if !args.no_header {
		 write!(handle, "taxon_id,taxon_name,taxon_rank");
		 if args.all_ranks {
			 for rank in rank::Rank::ranks() {
				 let rank_name = rank.to_string().replace(" ", "_");
				 write!(handle, ",{}_id,{}_name", rank_name, rank_name);
			 }
		 }
		 write!(handle, "\n");
	 }

	 for line in stdin.lock().lines() {
		  let line = line?;
		  if line.starts_with('>') {
				writeln!(handle, "{}", line)?;
		  } else {
				let id = line.parse::<taxon::TaxonId>()?;
				// Map to root if id not found
				let taxon = by_id.get_or_unknown(id)?;
				write!(handle, "{},{},{}", taxon.id, taxon.name, taxon.rank);
				if args.all_ranks {
					let lineage = by_id.lineage(id)?;
					for rank in rank::Rank::ranks() {
						if let Some(l_taxon) = &lineage[rank] {
							write!(handle, ",{},{}", l_taxon.id, l_taxon.name);
						} else {
							write!(handle, ",,");
						}
					}
				}
				write!(handle, "\n");
		  }
	 }
	 Ok(())
}

fn snaprank(args: args::SnapRank) -> Result<()> {
	let taxons = taxon::read_taxa_file(&args.taxon_file)?;
	if args.rank == rank::Rank::NoRank {
		return Err(ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
	}

	// Parsing the taxons
	let tree	 = taxon::TaxonTree::new(&taxons);
	let by_id	= taxon::TaxonList::new(taxons);
	let snapping = tree.filter_ancestors(|tid|
		by_id.get(tid).map(|t| t.rank == args.rank).unwrap_or(false)
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

fn jsontree(args: args::JsonTree) -> Result<()> {
	let taxons = taxon::read_taxa_file(args.taxon_file)?;

	// Parsing the taxons
	let tree	 = taxon::TaxonTree::new(&taxons);
	let by_id	= taxon::TaxonList::new(taxons);
	let snapping = tree.snapping(&by_id, args.ranked_only);

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

fn seedextend(args: args::SeedExtend) -> Result<()> {
	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

	let by_id = if let Some(ref tf) = args.ranked {
		let taxa = taxon::read_taxa_file(tf)?;
		Some(taxon::TaxonList::new_with_unknown(taxa, true))
	} else { None };

	for record in fasta::Reader::new(io::stdin(), false).records() {
		let record = record?;
		let mut taxons = record.sequence.iter()
		                       .map(|s| s.parse::<TaxonId>())
		                       .collect::<std::result::Result<Vec<TaxonId>,_>>()?;
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
				if same_max >= args.min_seed_size { seeds.push((start, end - same_tid)) }
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
		if same_max >= args.min_seed_size {
			if last_tid == 0 { end -= same_tid }
			seeds.push((start, end))
		}

		if let Some(ref by_id) = by_id {
			seeds = seeds.into_iter()
			             .max_by_key(|(s, e)| taxons.iter().skip(*s).take(e - s)
			                                        .map(|t| by_id.score(*t).unwrap_or(args.penalty))
			                                        .sum::<usize>())
			             .into_iter()
			             .collect::<Vec<(usize, usize)>>();
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

fn report(args: args::Report) -> Result<()> {
	let taxons = taxon::read_taxa_file(&args.taxon_file)?;
	if args.rank == rank::Rank::NoRank {
		return Err(ErrorKind::InvalidInvocation("Snap to an actual rank.".into()).into());
	}

	// Parsing the taxons
	let tree	 = taxon::TaxonTree::new(&taxons);
	let by_id	= taxon::TaxonList::new(taxons);
	let snapping = tree.filter_ancestors(|tid|
		by_id.get(tid).map(|t| t.rank == args.rank).unwrap_or(false)
	);

	// Read and count taxon ranks
	let mut counts = HashMap::new();
	let stdin = io::stdin();
	for line in stdin.lock().lines() {
		let taxon = line?.parse::<taxon::TaxonId>()?;
		*counts.entry(snapping[taxon].unwrap_or(0)).or_insert(0) += 1;
	}

	let mut counts = counts.into_iter()
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

fn bestof(args: args::BestOf) -> Result<()> {
	let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

	// Combine frames and process them
	let mut chunk = Vec::with_capacity(args.frames);
	for record in fasta::Reader::new(io::stdin(), false).records() {
		let record = record?;
		if chunk.len() < args.frames - 1 {
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
