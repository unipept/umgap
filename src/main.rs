use clap::Parser;

use umgap::commands;
use umgap::errors::Result;

fn main() {
    if let Err(err) = run() {
        // Mirrors what error_chain's quick_main! used to print.
        eprintln!("Error: {}", err);
        let mut source = std::error::Error::source(&err);
        while let Some(cause) = source {
            eprintln!("Caused by: {}", cause);
            source = cause.source();
        }
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    match Opt::parse() {
        Opt::BestOf(args) => commands::bestof::bestof(args),
        Opt::BuildIndex(args) => commands::buildindex::buildindex(args),
        Opt::FastqToFasta(args) => commands::fastq2fasta::fastq2fasta(args),
        Opt::Filter(args) => commands::filter::filter(args),
        Opt::JoinKmers(args) => commands::joinkmers::joinkmers(args),
        Opt::PeptToLca(args) => commands::pept2lca::pept2lca(args),
        Opt::PrintIndex(args) => commands::printindex::printindex(args),
        Opt::ProtToKmer(args) => commands::prot2kmer::prot2kmer(args),
        #[cfg(target_family = "unix")]
        Opt::ProtToKmerToLca(args) => commands::prot2kmer2lca::prot2kmer2lca(args),
        Opt::ProtToTryp(args) => commands::prot2tryp::prot2tryp(args),
        Opt::ProtToTrypToLca(args) => commands::prot2tryp2lca::prot2tryp2lca(args),
        Opt::SeedExtend(args) => commands::seedextend::seedextend(args),
        Opt::SnapTaxon(args) => commands::snaptaxon::snaptaxon(args),
        Opt::SplitKmers(args) => commands::splitkmers::splitkmers(args),
        Opt::TaxaToAgg(args) => commands::taxa2agg::taxa2agg(args),
        Opt::TaxaToFreq(args) => commands::taxa2freq::taxa2freq(args),
        Opt::TaxaToTree(args) => commands::taxa2tree::taxa2tree(args),
        Opt::Taxonomy(args) => commands::taxonomy::taxonomy(args),
        Opt::Translate(args) => commands::translate::translate(args),
        Opt::Uniq(args) => commands::uniq::uniq(args),
    }
}

/// UMGAP is a collection of tools to be used in metagenomics analysis pipelines. Use the
/// `umgap-analyse.sh` script for some prebuild pipelines.
///
/// Throughout this documentation, the term peptides is used for both tryptic peptides and k-mers.
/// The term taxon ID refers to an identifier of a NCBI taxonomy (which should be the same version
/// in the whole pipeline).
#[rustfmt::skip]
#[derive(Debug, Parser)]
#[command(version, propagate_version = true)]
pub enum Opt {
    #[command(name = "bestof")] BestOf(commands::bestof::BestOf),
    #[command(name = "buildindex")] BuildIndex(commands::buildindex::BuildIndex),
    #[command(name = "fastq2fasta")] FastqToFasta(commands::fastq2fasta::FastqToFasta),
    #[command(name = "filter")] Filter(commands::filter::Filter),
    #[command(name = "joinkmers")] JoinKmers(commands::joinkmers::JoinKmers),
    #[command(name = "pept2lca")] PeptToLca(commands::pept2lca::PeptToLca),
    #[command(name = "printindex")] PrintIndex(commands::printindex::PrintIndex),
    #[cfg(target_family = "unix")] #[command(name = "prot2kmer2lca")] ProtToKmerToLca(commands::prot2kmer2lca::ProtToKmerToLca),
    #[command(name = "prot2kmer")] ProtToKmer(commands::prot2kmer::ProtToKmer),
    #[command(name = "prot2tryp2lca")] ProtToTrypToLca(commands::prot2tryp2lca::ProtToTrypToLca),
    #[command(name = "prot2tryp")] ProtToTryp(commands::prot2tryp::ProtToTryp),
    #[command(name = "seedextend")] SeedExtend(commands::seedextend::SeedExtend),
    #[command(name = "snaptaxon")] SnapTaxon(commands::snaptaxon::SnapTaxon),
    #[command(name = "splitkmers")] SplitKmers(commands::splitkmers::SplitKmers),
    #[command(name = "taxa2agg")] TaxaToAgg(commands::taxa2agg::TaxaToAgg),
    #[command(name = "taxa2freq")] TaxaToFreq(commands::taxa2freq::TaxaToFreq),
    #[command(name = "taxa2tree")] TaxaToTree(commands::taxa2tree::TaxaToTree),
    #[command(name = "taxonomy")] Taxonomy(commands::taxonomy::Taxonomy),
    #[command(name = "translate")] Translate(commands::translate::Translate),
    #[command(name = "uniq")] Uniq(commands::uniq::Uniq),
}
