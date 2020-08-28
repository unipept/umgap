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
        Opt::ProtToPept(args) => commands::prot2pept::prot2pept(args),
        Opt::ProtToTrypToLca(args) => commands::prot2tryp2lca::prot2tryp2lca(args),
        Opt::Report(args) => commands::report::report(args),
        Opt::SeedExtend(args) => commands::seedextend::seedextend(args),
        Opt::SnapTaxon(args) => commands::snaptaxon::snaptaxon(args),
        Opt::SplitKmers(args) => commands::splitkmers::splitkmers(args),
        Opt::TaxaToAgg(args) => commands::taxa2agg::taxa2agg(args),
        Opt::Taxonomy(args) => commands::taxonomy::taxonomy(args),
        Opt::Translate(args) => commands::translate::translate(args),
        Opt::Uniq(args) => commands::uniq::uniq(args),
        Opt::Visualize(args) => commands::visualize::visualize(args),
    }
});

/// UMGAP is a collection of tools to be used in a metagenomics analysis pipeline.
///
/// Call `umgap <command> -h` for additional help per command.
///
/// Use the `umgap-analyse.sh` script for some prebuild pipelines.
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
    #[structopt(name = "prot2kmer")] ProtToKmer(commands::prot2kmer::ProtToKmer),
    #[structopt(name = "prot2kmer2lca")] ProtToKmerToLca(commands::prot2kmer2lca::ProtToKmerToLca),
    #[structopt(name = "prot2pept")] ProtToPept(commands::prot2pept::ProtToPept),
    #[structopt(name = "prot2tryp2lca")] ProtToTrypToLca(commands::prot2tryp2lca::ProtToTrypToLca),
    #[structopt(name = "report")] Report(commands::report::Report),
    #[structopt(name = "seedextend")] SeedExtend(commands::seedextend::SeedExtend),
    #[structopt(name = "snaptaxon")] SnapTaxon(commands::snaptaxon::SnapTaxon),
    #[structopt(name = "splitkmers")] SplitKmers(commands::splitkmers::SplitKmers),
    #[structopt(name = "taxa2agg")] TaxaToAgg(commands::taxa2agg::TaxaToAgg),
    #[structopt(name = "taxonomy")] Taxonomy(commands::taxonomy::Taxonomy),
    #[structopt(name = "translate")] Translate(commands::translate::Translate),
    #[structopt(name = "uniq")] Uniq(commands::uniq::Uniq),
    #[structopt(name = "visualize")] Visualize(commands::visualize::Visualize),
}
