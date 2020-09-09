use error_chain::quick_main;

use structopt::StructOpt;

use umgap::commands;
use umgap::errors::Result;

quick_main!(|| -> Result<()> {
    match Opt::from_args() {
        Opt::BestOf(args) => commands::bestof::bestof(args),
        Opt::BuildIndex(args) => commands::buildindex::buildindex(args),
        Opt::FastqToFasta(args) => commands::fastq2fasta::fastq2fasta(args),
        Opt::Filter(args) => commands::filter::filter(args),
        Opt::JoinKmers(args) => commands::joinkmers::joinkmers(args),
        Opt::PeptToLca(args) => commands::pept2lca::pept2lca(args),
        Opt::PrintIndex(args) => commands::printindex::printindex(args),
        Opt::ProtToKmer(args) => commands::prot2kmer::prot2kmer(args),
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
});

/// UMGAP is a collection of tools to be used in metagenomics analysis pipelines. Use the
/// `umgap-analyse.sh` script for some prebuild pipelines.
///
/// Throughout this documentation, the term taxon ID refers to an identifier of the NCBI taxonomy.
#[cfg_attr(rustfmt, rustfmt_skip)]
#[derive(Debug, StructOpt)]
pub enum Opt {
    #[structopt(name = "bestof")] BestOf(commands::bestof::BestOf),
    #[structopt(name = "buildindex")] BuildIndex(commands::buildindex::BuildIndex),
    #[structopt(name = "fastq2fasta")] FastqToFasta(commands::fastq2fasta::FastqToFasta),
    #[structopt(name = "filter")] Filter(commands::filter::Filter),
    #[structopt(name = "joinkmers")] JoinKmers(commands::joinkmers::JoinKmers),
    #[structopt(name = "pept2lca")] PeptToLca(commands::pept2lca::PeptToLca),
    #[structopt(name = "printindex")] PrintIndex(commands::printindex::PrintIndex),
    #[structopt(name = "prot2kmer2lca")] ProtToKmerToLca(commands::prot2kmer2lca::ProtToKmerToLca),
    #[structopt(name = "prot2kmer")] ProtToKmer(commands::prot2kmer::ProtToKmer),
    #[structopt(name = "prot2tryp2lca")] ProtToTrypToLca(commands::prot2tryp2lca::ProtToTrypToLca),
    #[structopt(name = "prot2tryp")] ProtToTryp(commands::prot2tryp::ProtToTryp),
    #[structopt(name = "seedextend")] SeedExtend(commands::seedextend::SeedExtend),
    #[structopt(name = "snaptaxon")] SnapTaxon(commands::snaptaxon::SnapTaxon),
    #[structopt(name = "splitkmers")] SplitKmers(commands::splitkmers::SplitKmers),
    #[structopt(name = "taxa2agg")] TaxaToAgg(commands::taxa2agg::TaxaToAgg),
    #[structopt(name = "taxa2freq")] TaxaToFreq(commands::taxa2freq::TaxaToFreq),
    #[structopt(name = "taxa2tree")] TaxaToTree(commands::taxa2tree::TaxaToTree),
    #[structopt(name = "taxonomy")] Taxonomy(commands::taxonomy::Taxonomy),
    #[structopt(name = "translate")] Translate(commands::translate::Translate),
    #[structopt(name = "uniq")] Uniq(commands::uniq::Uniq),
}
