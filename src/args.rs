//! Argument parsing for the UMGAP

use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;

use rank::Rank;
use taxon::TaxonId;

/// A reading frame
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy)]
pub enum Frame {
	Forward1,
	Forward2,
	Forward3,
	Reverse1,
	Reverse2,
	Reverse3,
}

static FRAMES: &[&str] = &["1", "2", "3", "1R", "2R", "3R"];
impl Frame {
	fn variants() -> &'static [&'static str] {
		FRAMES
	}
}

impl FromStr for Frame {
	type Err = Error;

	fn from_str(s: &str) -> Result<Self> {
		match s {
			"1" => Ok(Frame::Forward1),
			"2" => Ok(Frame::Forward2),
			"3" => Ok(Frame::Forward3),
			"1R" => Ok(Frame::Reverse1),
			"2R" => Ok(Frame::Reverse2),
			"3R" => Ok(Frame::Reverse3),
			_ => Err(ErrorKind::ParseFrameError(s.to_string()).into()),
		}
	}
}

impl fmt::Display for Frame {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "{}", match *self {
			Frame::Forward1 => "1",
			Frame::Forward2 => "2",
			Frame::Forward3 => "3",
			Frame::Reverse1 => "1R",
			Frame::Reverse2 => "2R",
			Frame::Reverse3 => "3R",
		})
	}
}

/// An aggregation method
#[allow(missing_docs)]
#[derive(Debug)]
pub enum Method {
	Tree,
	RangeMinimumQuery,
}

static METHODS: &[&str] = &["tree", "rmq"];
impl Method {
	fn variants() -> &'static [&'static str] {
		METHODS
	}
}

impl FromStr for Method {
	type Err = Error;

	fn from_str(s: &str) -> Result<Self> {
		match s {
			"tree" => Ok(Method::Tree),
			"rmq" => Ok(Method::RangeMinimumQuery),
			_ => Err(ErrorKind::ParseMethodError(s.to_string()).into()),
		}
	}
}

/// An aggregation strategy
#[allow(missing_docs)]
#[derive(Debug)]
pub enum Strategy {
	LowestCommonAncestor,
	Hybrid,
	MaximumRootToLeafPath,
}

static STRATEGIES: &[&str] = &["lca*", "hybrid", "mrtl"];
impl Strategy {
	fn variants() -> &'static [&'static str] {
		STRATEGIES
	}
}

impl FromStr for Strategy {
	type Err = Error;

	fn from_str(s: &str) -> Result<Self> {
		match s {
			"lca*" => Ok(Strategy::LowestCommonAncestor),
			"hybrid" => Ok(Strategy::Hybrid),
			"mrtl" => Ok(Strategy::MaximumRootToLeafPath),
			_ => Err(ErrorKind::ParseStrategyError(s.to_string()).into()),
		}
	}
}

/// The Options enum for UMGAP arguments
#[derive(Debug, StructOpt)]
#[allow(missing_docs)]
pub enum Opt {
	/// Translates DNA into Amino Acid Sequences.
	#[structopt(name = "translate")]
	Translate(Translate),

	/// Looks up each line of input in a given FST index and outputs
	/// the result. Lines starting with a '>' are copied. Lines for
	/// which no mapping is found are ignored.
	#[structopt(name = "pept2lca")]
	PeptToLca(PeptToLca),

	/// Reads all the records in a specified FASTA file and queries the
	/// k-mers in an FST for the LCA's.
	#[structopt(name = "prot2kmer2lca")]
	ProtToKmerToLca(ProtToKmerToLca),

	/// Reads all the records in a specified FASTA file and queries the
	/// tryptic peptides in an FST for the LCA's.
	#[structopt(name = "prot2tryp2lca")]
	ProtToTrypToLca(ProtToTrypToLca),

	/// Aggregates taxa to a single taxon.
	#[structopt(name = "taxa2agg")]
	TaxaToAgg(TaxaToAgg),

	/// Splits each protein sequence in a FASTA format into a list of (tryptic) peptides.
	#[structopt(name = "prot2pept")]
	ProtToPept(ProtToPept),

	/// Pick the frame with the most none-root hits.
	#[structopt(name = "bestof")]
	BestOf(BestOf),

	/// Count and report on a list of taxon ids.
	#[structopt(name = "report")]
	Report(Report),

	/// Seed and extend.
	#[structopt(name = "seedextend")]
	SeedExtend(SeedExtend),

	/// Aggregates taxa to a JSON tree for usage in the unipept visualisations.
	#[structopt(name = "jsontree")]
	JsonTree(JsonTree),

	/// Show taxonomy info
	#[structopt(name = "taxonomy")]
	Taxonomy(Taxonomy),

	/// Snap taxa to a specified rank or one of the specified taxa.
	#[structopt(name = "snaptaxon")]
	SnapTaxon(SnapTaxon),

	/// Interleaves a number of FASTQ files into a single FASTA output.
	#[structopt(name = "fastq2fasta")]
	FastqToFasta(FastqToFasta),

	/// Filter peptides in a FASTA format based on specific criteria.
	#[structopt(name = "filter")]
	Filter(Filter),

	/// Concatenates the data strings of all consecutive FASTA entries with the same header.
	#[structopt(name = "uniq")]
	Uniq(Uniq),

	/// Splits each protein sequence in a FASTA format into a list of kmers.
	#[structopt(name = "prot2kmer")]
	ProtToKmer(ProtToKmer),

	/// Print the values in an FST index to stdout.
	#[structopt(name = "printindex")]
	PrintIndex(PrintIndex),

	/// Write an FST index of stdin on stdout.
	#[structopt(name = "buildindex")]
	BuildIndex,

	/// Count the amount of FASTA records from stdin
	#[structopt(name = "countrecords")]
	CountRecords,

	/// Read FASTA-records from stdin and create a frequency table from their
	/// sequences.
	#[structopt(name = "counts")]
	Counts,

	/// Look up FASTA-sequences in an index file
	#[structopt(name = "lookup")]
	Lookup(Lookup),

	/// Split FASTA-sequences into kmers and look them up in an index file
	#[structopt(name = "kmerlookup")]
	KmerLookup(KmerLookup),

	/// Aggregate results per record by picking the sequence with the most
	/// occurences.
	#[structopt(name = "aggregate")]
	Aggregate(Aggregate),
}

/// Translates DNA into Amino Acid Sequences.
#[derive(Debug, StructOpt)]
pub struct Translate {
	/// Replace each start-codon with methionine
	#[structopt(short = "m", long = "methionine")]
	pub methionine: bool,

	/// Read and output all six frames
	#[structopt(short = "a", long = "all-frames", conflicts_with = "frame")]
	pub all_frames: bool,

	/// Adds a reading frame (1, 2, 3, 1R, 2R or 3R)
	#[structopt(
	            short = "f",
	            long = "frame",
	            raw(possible_values = "&Frame::variants()")
	)]
	pub frames: Vec<Frame>,

	/// Append a bar (|) and the name of the frame to the fasta header
	#[structopt(short = "n", long = "append-name")]
	pub append_name: bool,

	/// Translation table to use
	#[structopt(short = "t", long = "table", default_value = "1")]
	pub table: String,

	/// Print the selected table and exit
	#[structopt(short = "s", long = "show-table")]
	pub show_table: bool,
}

/// Looks up each line of input in a given FST index and outputs
/// the result. Lines starting with a '>' are copied. Lines for
/// which no mapping is found are ignored.
#[derive(Debug, StructOpt)]
pub struct PeptToLca {
	/// Map unknown sequences to 0 instead of ignoring them
	#[structopt(short = "o", long = "one-on-one")]
	pub one_on_one: bool,

	/// An FST to query
	#[structopt(parse(from_os_str))]
	pub fst_file: PathBuf,

	/// Load FST in memory instead of mmap'ing the file contents. This makes
	/// querying significantly faster, but requires some time to load the FST
	/// into memory.
	#[structopt(short = "m", long = "in-memory")]
	pub fst_in_memory: bool,

	/// How much records to group together in one chunk, bigger chunks decrease
	/// the overhead caused by multithreading. Because the output order is not
	/// necessarily the same as the input order, having a chunk size which is
	/// a multiple of 12 (all 6 translations multiplied by the two paired-end
	/// reads) will keep sequences of the same reads together.
	#[structopt(short = "c", long = "chunksize", default_value = "240")]
	pub chunk_size: usize,

	/// How much threads can be used for parallelization. Default is to let the
	/// underlying library [rayon](https://docs.rs/rayon/latest/rayon/) select
	/// the thread count automatically based on the number of logical CPUs.
	#[structopt(short = "j", long = "thread-count")]
	pub thread_count: Option<usize>,
}

/// Reads all the records in a specified FASTA file and queries the
/// k-mers in an FST for the LCA's.
#[derive(Debug, StructOpt)]
pub struct ProtToKmerToLca {
	/// The length of the k-mers in the FST
	#[structopt(short = "k", long = "length", default_value = "9")]
	pub length: usize,

	/// Map unknown sequences to 0 instead of ignoring them
	#[structopt(short = "o", long = "one-on-one")]
	pub one_on_one: bool,

	/// An FST to query
	#[structopt(parse(from_os_str))]
	pub fst_file: PathBuf,

	/// Instead of reading from stdin and writing to stdout, create a Unix
	/// socket to communicate with. This socket can be communicated with using
	/// socat for example: (` socat - UNIX-CONNECT:<socket>`).
	/// This is especially useful in combination with the `--in-memory` flag:
	/// you only have to load the FST in memory once, after which you can query
	/// it without having the loading time overhead each time.
	#[structopt(parse(from_os_str), short = "s", long = "socket")]
	pub socket: Option<PathBuf>,

	/// Load FST in memory instead of mmap'ing the file contents. This makes
	/// querying significantly faster, but requires some time to load the FST
	/// into memory.
	#[structopt(short = "m", long = "in-memory")]
	pub fst_in_memory: bool,

	/// How much records to group together in one chunk, bigger chunks decrease
	/// the overhead caused by multithreading. Because the output order is not
	/// necessarily the same as the input order, having a chunk size which is
	/// a multiple of 12 (all 6 translations multiplied by the two paired-end
	/// reads) will keep sequences of the same reads together.
	#[structopt(short = "c", long = "chunksize", default_value = "240")]
	pub chunk_size: usize,

	/// How much threads can be used for parallelization. Default is to let the
	/// underlying library [rayon](https://docs.rs/rayon/latest/rayon/) select
	/// the thread count automatically based on the number of logical CPUs.
	#[structopt(short = "j", long = "thread-count")]
	pub thread_count: Option<usize>,
}

/// Reads all the records in a specified FASTA file and queries the
/// tryptic peptides in an FST for the LCA's.
#[derive(Debug, StructOpt)]
pub struct ProtToTrypToLca {
	/// Map unknown sequences to 0 instead of ignoring them
	#[structopt(short = "o", long = "one-on-one")]
	pub one_on_one: bool,

	/// An FST to query
	#[structopt(parse(from_os_str))]
	pub fst_file: PathBuf,

	/// Load FST in memory instead of mmap'ing the file contents. This makes
	/// querying significantly faster, but requires some time to load the FST
	/// into memory.
	#[structopt(short = "m", long = "in-memory")]
	pub fst_in_memory: bool,

	/// How much reads to group together in one chunk, bigger chunks decrease
	/// the overhead caused by multithreading. Because the output order is not
	/// necessarily the same as the input order, having a chunk size which is
	/// a multiple of 12 (all 6 translations multiplied by the two paired-end
	/// reads) will keep sequences of the same reads together.
	#[structopt(short = "c", long = "chunksize", default_value = "240")]
	pub chunk_size: usize,

	/// The cleavage-pattern (regex), i.e. the pattern after which
	/// the next peptide will be cleaved for tryptic peptides)
	#[structopt(short = "p", long = "pattern", default_value = "([KR])([^P])")]
	pub pattern: String,

	/// Minimum length of tryptic peptides to be mapped
	#[structopt(short = "l", long = "minlen", default_value = "5")]
	pub min_length: usize,

	/// Maximum length of tryptic peptides to be mapped
	#[structopt(short = "L", long = "maxlen", default_value = "50")]
	pub max_length: usize,

	/// The letters that a tryptic peptide must contain
	#[structopt(short = "c", long = "contains", default_value = "")]
	pub contains: String,

	/// The letters that a tryptic peptide mustn't contain
	#[structopt(short = "l", long = "lacks", default_value = "")]
	pub lacks: String,
}

/// Aggregates taxa to a single taxon.
#[derive(Debug, StructOpt)]
pub struct TaxaToAgg {
	/// Each taxon is followed by a score between 0 and 1
	#[structopt(short = "s", long = "scored")]
	pub scored: bool,

	/// Restrict to taxa with a taxonomic rank
	#[structopt(short = "r", long = "ranked")]
	pub ranked_only: bool,

	/// The method to use for aggregation
	#[structopt(
	            short = "m",
	            long = "method",
	            default_value = "tree",
	            raw(possible_values = "&Method::variants()")
	)]
	pub method: Method,

	/// The method to use for aggregation
	#[structopt(
	            short = "a",
	            long = "aggregate",
	            default_value = "hybrid",
	            raw(possible_values = "&Strategy::variants()")
	)]
	pub strategy: Strategy,

	/// The factor for the hybrid aggregation, from 0.0 (MRTL) to
	/// 1.0 (LCA*)
	#[structopt(short = "f", long = "factor", default_value = "0.25")]
	pub factor: f32,

	/// The smallest input frequency for a taxon to be included in
	/// the aggregation
	#[structopt(short = "l", long = "lower-bound", default_value = "0.0")]
	pub lower_bound: f32,

	/// The NCBI taxonomy tsv-file
	#[structopt(parse(from_os_str))]
	pub taxon_file: PathBuf,
}

/// Splits each protein sequence in a FASTA format into a list of (tryptic) peptides.
#[derive(Debug, StructOpt)]
pub struct ProtToPept {
	/// The cleavage-pattern (regex), i.e. the pattern after which
	/// the next peptide will be cleaved for tryptic peptides)
	#[structopt(short = "p", long = "pattern", default_value = "([KR])([^P])")]
	pub pattern: String,
}

/// Splits each protein sequence in a FASTA format into a list of kmers.
#[derive(Debug, StructOpt)]
pub struct ProtToKmer {
	/// The K in K-mers
	#[structopt(short = "k", long = "length", default_value = "9")]
	pub length: usize,
}

/// Concatenates the data strings of all consecutive FASTA entries with the same header.
#[derive(Debug, StructOpt)]
pub struct Uniq {
	/// Separator between output items
	#[structopt(short = "s", long = "separator", default_value = "\n")]
	pub separator: String,

	/// Wrap the output sequences
	#[structopt(short = "w", long = "wrap")]
	pub wrap: bool,
}

/// Filter peptides in a FASTA format based on specific criteria.
#[derive(Debug, StructOpt)]
pub struct Filter {
	/// Minimum length required
	#[structopt(short = "m", long = "minlen", default_value = "5")]
	pub min_length: usize,

	/// Maximum length allowed
	#[structopt(short = "M", long = "maxlen", default_value = "50")]
	pub max_length: usize,

	/// The letters that a sequence must contain
	#[structopt(short = "c", long = "contains", default_value = "")]
	pub contains: String,

	/// The letters that a sequence mustn't contain
	#[structopt(short = "l", long = "lacks", default_value = "")]
	pub lacks: String,
}

/// Interleaves a number of FASTQ files into a single FASTA output.
#[derive(Debug, StructOpt)]
pub struct FastqToFasta {
	/// The input files
	#[structopt(parse(from_os_str))]
	pub input: Vec<PathBuf>,
}

/// Snap taxa to a specified rank or one of the specified taxa.
#[derive(Debug, StructOpt)]
pub struct SnapTaxon {
	/// The rank to snap towards.
	#[structopt(short = "r", long = "rank")]
	pub rank: Option<Rank>,

	/// A taxon to snap towards (allow multiple times).
	#[structopt(short = "t", long = "taxons")]
	pub taxons: Vec<TaxonId>,

	/// Include invalidated taxa
	#[structopt(short = "i", long = "invalid")]
	pub invalid: bool,

	/// The NCBI taxonomy tsv-file
	#[structopt(parse(from_os_str))]
	pub taxon_file: PathBuf,
}

/// Show the id, name, and rank of the given taxon ids in a CSV format
#[derive(Debug, StructOpt)]
pub struct Taxonomy {
	/// The NCBI taxonomy tsv-file
	#[structopt(parse(from_os_str))]
	pub taxon_file: PathBuf,

	/// Show the full lineage of a taxon. Ranks below the given taxon
	/// whill be empty.
	#[structopt(short = "a", long = "all")]
	pub all_ranks: bool,

	/// Do not output the CSV header
	#[structopt(short = "H", long = "no-header")]
	pub no_header: bool,
}


/// Aggregates taxa to a JSON tree for usage in the unipept visualisations.
#[derive(Debug, StructOpt)]
pub struct JsonTree {
	/// Exclude taxons with unnamed rank
	#[structopt(short = "r", long = "ranked")]
	pub ranked_only: bool,

	/// The minimum frequency to be reported
	#[structopt(short = "f", long = "frequency", default_value = "1")]
	pub min_frequency: usize,

	/// The NCBI taxonomy tsv-file
	#[structopt(parse(from_os_str))]
	pub taxon_file: PathBuf,
}

/// Seed and extend.
#[derive(Debug, StructOpt)]
pub struct SeedExtend {
	/// The minimum length of equal taxa to count as seed
	#[structopt(short = "s", long = "min-seed-size", default_value = "2")]
	pub min_seed_size: usize,

	/// The maximum length of a gap between seeds in an extension
	#[structopt(short = "g", long = "max-gap-size", default_value = "0")]
	pub max_gap_size: usize,

	/// Use taxon ranks in given NCBI taxonomy tsv-file to pick extended seed with highest score
	#[structopt(short = "r", long = "ranked", parse(from_os_str))]
	pub ranked: Option<PathBuf>,

	/// The score penalty for gaps in extended seeds
	#[structopt(short = "p", long = "penalty", default_value = "5")]
	pub penalty: usize,
}

/// Count and report on a list of taxon ids.
#[derive(Debug, StructOpt)]
pub struct Report {
	/// The rank to show
	#[structopt(short = "r", long = "rank", default_value = "species")]
	pub rank: Rank,

	/// The minimum frequency to be reported
	#[structopt(short = "f", long = "frequency", default_value = "1")]
	pub min_frequency: usize,

	/// The NCBI taxonomy tsv-file
	#[structopt(parse(from_os_str))]
	pub taxon_file: PathBuf,
}

/// Pick the frame with the most none-root hits.
#[derive(Debug, StructOpt)]
pub struct BestOf {
	/// The number of frames of which to pick the best
	#[structopt(short = "f", long = "frames", default_value = "6")]
	pub frames: usize,
}

/// Print the values in an FST index to stdout.
#[derive(Debug, StructOpt)]
pub struct PrintIndex {
	/// An FST to query
	#[structopt(parse(from_os_str))]
	pub fst_file: PathBuf,
}

/// Use the sequences of FASTA-records in stdin to look them up
/// in the index file.
#[derive(Debug, StructOpt)]
pub struct Lookup {
	/// Index file to query. The format of this file is <key><delimiter><value>.
	/// <key> kan be any string, but each key has to be unique.
	/// <delimiter> can be set with -d (default is tab).
	/// <value> can be any string and may contain <delimiter>.
	#[structopt(parse(from_os_str))]
	pub index_file: PathBuf,

	/// Which delimiter the index file uses to separate keys from values
	/// (default: tab).
	#[structopt(short = "d", long = "delimiter", default_value = "\t")]
	pub delimiter: String,

	/// How much records to group together in one chunk, bigger chunks decrease
	/// the overhead caused by multithreading. Because the output order is not
	/// necessarily the same as the input order, having a chunk size which is
	/// a multiple of 12 (all 6 translations multiplied by the two paired-end
	/// reads) will keep sequences of the same reads together.
	#[structopt(short = "c", long = "chunksize", default_value = "240")]
	pub chunk_size: usize,

	/// How much threads can be used for parallelization. Default is to let the
	/// underlying library [rayon](https://docs.rs/rayon/latest/rayon/) select
	/// the thread count automatically based on the number of logical CPUs.
	#[structopt(short = "j", long = "thread-count")]
	pub thread_count: Option<usize>,

	/// Map unknown keys to the given default value instead of ignoring them
	#[structopt(short = "f", long = "default")]
	pub default: Option<String>,

	/// Skip first line of the index file
	#[structopt(short = "h", long = "header")]
	pub has_header: bool,
}

/// Split the sequences of FASTA-records in stdin into kmers and look them up
/// in the index file
#[derive(Debug, StructOpt)]
pub struct KmerLookup {
	/// Index file to query. The format of this file is <key><delimiter><value>.
	/// <key> kan be any string with a length of <kmer_length>, but each key
	/// has to be unique.
	/// <delimiter> can be set with -d (default is tab).
	/// <value> can be any string and may contain <delimiter>.
	#[structopt(parse(from_os_str))]
	pub index_file: PathBuf,

	/// Which delimiter the index file uses to separate keys from values
	/// (default: tab).
	#[structopt(short = "d", long = "delimiter", default_value = "\t")]
	pub delimiter: String,

	/// How much records to group together in one chunk, bigger chunks decrease
	/// the overhead caused by multithreading. Because the output order is not
	/// necessarily the same as the input order, having a chunk size which is
	/// a multiple of 12 (all 6 translations multiplied by the two paired-end
	/// reads) will keep sequences of the same reads together.
	#[structopt(short = "c", long = "chunksize", default_value = "240")]
	pub chunk_size: usize,

	/// How much threads can be used for parallelization. Default is to let the
	/// underlying library [rayon](https://docs.rs/rayon/latest/rayon/) select
	/// the thread count automatically based on the number of logical CPUs.
	#[structopt(short = "j", long = "thread-count")]
	pub thread_count: Option<usize>,

	/// Map unknown keys to the given default value instead of ignoring them
	#[structopt(short = "f", long = "default")]
	pub default: Option<String>,

	/// Skip first line of the index file
	#[structopt(short = "h", long = "header")]
	pub has_header: bool,

	/// The K in K-mers
	#[structopt(short = "k", long = "length", default_value = "9")]
	pub kmer_length: usize,

	/// Instead of reading from stdin and writing to stdout, create a Unix
	/// socket to communicate with. This socket can be communicated with using
	/// socat for example: (` socat - UNIX-CONNECT:<socket>`).
	#[structopt(parse(from_os_str), short = "s", long = "socket")]
	pub socket: Option<PathBuf>,
}

/// Aggregate results per record by picking the sequence with the most
/// occurences. When multiple sequences occur the most, one is arbitrarily (but
/// deterministically) selected.
#[derive(Debug, StructOpt)]
pub struct Aggregate {
	/// Minimum amount of occurences before a sequence is selected
	#[structopt(short = "m", long = "minimum", default_value = "0")]
	pub minimum: usize,
}


error_chain! {
	errors {
		/// Unparseable Frame
		ParseFrameError(frame: String) {
			description("Unparseable frame")
			display("Unparseable frame: {}", frame)
		}
		/// Unparseable Method
		ParseMethodError(method: String) {
			description("Unparseable method")
			display("Unparseable method: {}", method)
		}
		/// Unparseable Strategy
		ParseStrategyError(strategy: String) {
			description("Unparseable strategy")
			display("Unparseable strategy: {}", strategy)
		}
	}
}
