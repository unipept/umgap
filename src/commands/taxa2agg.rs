//! The `umgap taxa2agg` command.

use std::io;
use std::io::Write;
use std::path::PathBuf;
use std::str::FromStr;

use crate::agg;
use crate::errors;
use crate::io::fasta;
use crate::rmq;
use crate::taxon;
use crate::taxon::TaxonId;
use crate::tree;

/// The `umgap taxa2agg` command takes one or more lists of taxa as input and reduces each to a
/// single consensus taxon.
///
/// The input is given on *standard input*, in a FASTA format. Per FASTA header, there should be a
/// list of taxa, each represented by its identifier, seperated by newlines. The output is written
/// to *standard output*, also in a FASTA format. Each header in the output is accompanied by a
/// single taxon identifier, which is the consensus taxon resulting from aggregation of the given
/// list.
///
/// For operation, the command requires to be passed a small database created from the uniprot
/// taxonomy as argument.
///
///     $ cat input.fa
///     >header1
///     571525
///     571525
///     6920
///     6920
///     1
///     6920
///     $ umgap taxa2agg taxons.tsv < input.fa
///     >header1
///     571525
///
/// By default, the aggregation used is the maximum root-to-leaf path (MRTL) aggregation. A variant
/// of the lowest common ancestor (LCA\*) aggregation is also available via the `-a` and `-m`
/// options, as is a hybrid approach.
///
/// * `-m rmq -a mrtl` is the default aggregation strategy. It selects the taxon from the given list
///   which has the highest frequency of ancestors in the list (including its own frequency). A
///   range-minimum-query (RMQ) method is used.
///
/// * `-m tree -a lca\*` returns the taxon (possibly not from the list) of lowest rank without
///   contradicting taxa in the list. Non-contradicting taxa of a taxon are either itself, its
///   ancestors and its descendants. A tree-based method is used.
///
/// * `-m tree -a hybrid` mixes above two strategies, which results in a taxon which might have not
///   have the highest frequency of ancestors in the list, but would have less contradicting taxa.
///   Use the `-f` option to select a hybrid close to the MTRL (`-f 0.0`) or to the LCA (`-f 1.0`).
#[derive(Debug, StructOpt)]
pub struct TaxaToAgg {
    /// Each taxon is followed by a score between 0 and 1
    #[structopt(short = "s", long = "scored")]
    pub scored: bool,

    /// Let all taxa snap to taxa with a named rank (such as species) during calculations
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

    /// The strategy to use for aggregation
    #[structopt(
                short = "a",
                long = "aggregate",
                default_value = "hybrid",
                possible_values = &Strategy::variants()
    )]
    pub strategy: Strategy,

    /// The factor for the hybrid aggregation, from 0.0 (MRTL) to 1.0 (LCA*)
    #[structopt(short = "f", long = "factor", default_value = "0.25")]
    pub factor: f32,

    /// The smallest input frequency for a taxon to be included in the aggregation
    #[structopt(short = "l", long = "lower-bound", default_value = "0")]
    pub lower_bound: f32,

    /// The NCBI taxonomy tsv-file
    #[structopt(parse(from_os_str))]
    pub taxon_file: PathBuf,
}

/// Implements the taxa2agg command.
pub fn taxa2agg(args: TaxaToAgg) -> errors::Result<()> {
    // Parsing the Taxa file
    let taxons = taxon::read_taxa_file(args.taxon_file)?;

    // Parsing the taxons
    let tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = tree.snapping(&by_id, args.ranked_only);

    let aggregator: errors::Result<Box<dyn agg::Aggregator>> = match (args.method, args.strategy) {
        (Method::RangeMinimumQuery, Strategy::MaximumRootToLeafPath) => {
            Ok(Box::new(rmq::rtl::RTLCalculator::new(tree.root, &by_id)))
        }
        (Method::RangeMinimumQuery, Strategy::LowestCommonAncestor) => {
            Ok(Box::new(rmq::lca::LCACalculator::new(tree)))
        }
        (Method::RangeMinimumQuery, Strategy::Hybrid) => {
            writeln!(
                &mut io::stderr(),
                "Warning: this is a hybrid between LCA/MRTL, not LCA*/MRTL"
            )
            .unwrap();
            Ok(Box::new(rmq::mix::MixCalculator::new(tree, args.factor)))
        }
        (Method::Tree, Strategy::LowestCommonAncestor) => {
            Ok(Box::new(tree::lca::LCACalculator::new(tree.root, &by_id)))
        }
        (Method::Tree, Strategy::Hybrid) => Ok(Box::new(tree::mix::MixCalculator::new(
            tree.root,
            &by_id,
            args.factor,
        ))),
        (m, s) => Err(errors::ErrorKind::InvalidInvocation(format!(
            "{:?} and {:?} cannot be combined",
            m, s
        ))
        .into()),
    };
    let aggregator = aggregator?;

    fn with_score(pair: &String) -> errors::Result<(TaxonId, f32)> {
        let split = pair.split('=').collect::<Vec<_>>();
        if split.len() != 2 {
            Err("Taxon without score")?;
        }
        Ok((split[0].parse::<TaxonId>()?, split[1].parse::<f32>()?))
    }

    fn not_scored(tid: &String) -> errors::Result<(TaxonId, f32)> {
        Ok((tid.parse::<TaxonId>()?, 1.0))
    }

    let parser = if args.scored { with_score } else { not_scored };

    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    // Iterate over each read
    for record in fasta::Reader::new(io::stdin(), false).records() {
        // Parse the sequence of LCA's
        let record = record?;
        let taxons = record
            .sequence
            .iter()
            .map(parser)
            .collect::<errors::Result<Vec<(TaxonId, f32)>>>()?;

        // Create a frequency table of taxons for this read (taking into account the lower bound)
        let counts = agg::count(taxons.into_iter().filter(|&(tid, _)| tid != 0));
        let counts = agg::filter(counts, args.lower_bound);

        writer.write_record(fasta::Record {
            header: record.header,
            sequence: if counts.is_empty() {
                vec!["1".into()]
            } else {
                let aggregate = aggregator.aggregate(&counts)?;
                vec![snapping[aggregate].unwrap().to_string()]
            },
        })?;
    }
    Ok(())
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
