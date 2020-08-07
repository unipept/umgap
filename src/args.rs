//! Argument parsing for the UMGAP

use std::path::PathBuf;
use std::str::FromStr;

use crate::commands;
use crate::rank::Rank;
use crate::taxon::TaxonId;

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
    Translate(commands::translate::Translate),

    /// Looks up each line of input in a given FST index and outputs
    /// the result. Lines starting with a '>' are copied. Lines for
    /// which no mapping is found are ignored.
    #[structopt(name = "pept2lca")]
    PeptToLca(commands::pept2lca::PeptToLca),

    /// Reads all the records in a specified FASTA file and queries the
    /// k-mers in an FST for the LCA's.
    #[structopt(name = "prot2kmer2lca")]
    ProtToKmerToLca(commands::prot2kmer2lca::ProtToKmerToLca),

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

    /// Splits each protein sequence in a FASTA format into a list of kmers.
    #[structopt(name = "splitkmers")]
    SplitKmers(SplitKmers),

    /// Groups a CSV by equal first fields (Kmers) and aggregates the second fields (taxon ids).
    #[structopt(name = "joinkmers")]
    JoinKmers(JoinKmers),

    /// Write an FST index of stdin on stdout.
    #[structopt(name = "buildindex")]
    BuildIndex,

    /// Count the amount of FASTA records from stdin
    #[structopt(name = "countrecords")]
    CountRecords,

    /// Visualizes the given list of taxons using the Unipept API
    #[structopt(name = "visualize")]
    Visualize(Visualize),
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

    /// Keep only tryptic peptides containing these letters
    #[structopt(short = "k", long = "keep", default_value = "")]
    pub contains: String,

    /// Drop tryptic peptides containing these letters
    #[structopt(short = "d", long = "drop", default_value = "")]
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
	            possible_values = &Method::variants()
	)]
    pub method: Method,

    /// The method to use for aggregation
    #[structopt(
	            short = "a",
	            long = "aggregate",
	            default_value = "hybrid",
	            possible_values = &Strategy::variants()
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

/// Splits each taxon id + protein sequence pair in a CSV format into a list of kmers.
#[derive(Debug, StructOpt)]
pub struct SplitKmers {
    /// The K in K-mers
    #[structopt(short = "k", long = "length", default_value = "9")]
    pub length: usize,
    /// Print (K-1)-mer suffixes of the Kmers starting with this character
    #[structopt(short = "p", long = "prefix", default_value = "")]
    pub prefix: String,
}

/// Groups a CSV by equal first fields (Kmers) and aggregates the second fields (taxon ids).
#[derive(Debug, StructOpt)]
pub struct JoinKmers {
    /// The NCBI taxonomy tsv-file
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,
}

/// Visualizes the given list of taxons using the Unipept API
#[derive(Debug, StructOpt)]
pub struct Visualize {
    /// Host the result online and return the URL
    #[structopt(short = "u", long = "url")]
    pub url: bool,
}

error_chain! {
    errors {
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
