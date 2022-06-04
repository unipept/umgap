use error_chain::quick_main;

use structopt::StructOpt;

use umgap::cli;
use umgap::errors::Result;

quick_main!(|| -> Result<()> {
    match Opt::from_args() {
        Opt::BestOf(args) => cli::bestof::bestof(args),
        Opt::BuildIndex(args) => cli::buildindex::buildindex(args),
        Opt::FastqToFasta(args) => cli::fastq2fasta::fastq2fasta(args),
        Opt::Filter(args) => cli::filter::filter(args),
        Opt::JoinKmers(args) => cli::joinkmers::joinkmers(args),
        Opt::PeptToLca(args) => cli::pept2lca::pept2lca(args),
        Opt::PrintIndex(args) => cli::printindex::printindex(args),
        Opt::ProtToKmer(args) => cli::prot2kmer::prot2kmer(args),
        #[cfg(target_family = "unix")]
        Opt::ProtToKmerToLca(args) => cli::prot2kmer2lca::prot2kmer2lca(args),
        Opt::ProtToTryp(args) => cli::prot2tryp::prot2tryp(args),
        Opt::ProtToTrypToLca(args) => cli::prot2tryp2lca::prot2tryp2lca(args),
        Opt::SeedExtend(args) => cli::seedextend::seedextend(args),
        Opt::SnapTaxon(args) => cli::snaptaxon::snaptaxon(args),
        Opt::SplitKmers(args) => cli::splitkmers::splitkmers(args),
        Opt::TaxaToAgg(args) => cli::taxa2agg::taxa2agg(args),
        Opt::TaxaToFreq(args) => cli::taxa2freq::taxa2freq(args),
        Opt::TaxaToTree(args) => cli::taxa2tree::taxa2tree(args),
        Opt::Taxonomy(args) => cli::taxonomy::taxonomy(args),
        Opt::Translate(args) => cli::translate::translate(args),
        Opt::Uniq(args) => cli::uniq::uniq(args),
    }
});

/// UMGAP is a collection of tools to be used in metagenomics analysis pipelines. Use the
/// `umgap-analyse.sh` script for some prebuild pipelines.
///
/// Throughout this documentation, the term peptides is used for both tryptic peptides and k-mers.
/// The term taxon ID refers to an identifier of a NCBI taxonomy (which should be the same version
/// in the whole pipeline).
#[rustfmt::skip]
#[derive(Debug, StructOpt)]
pub enum Opt {
    #[structopt(name = "bestof")] BestOf(cli::bestof::BestOf),
    #[structopt(name = "buildindex")] BuildIndex(cli::buildindex::BuildIndex),
    #[structopt(name = "fastq2fasta")] FastqToFasta(cli::fastq2fasta::FastqToFasta),
    #[structopt(name = "filter")] Filter(cli::filter::Filter),
    #[structopt(name = "joinkmers")] JoinKmers(cli::joinkmers::JoinKmers),
    #[structopt(name = "pept2lca")] PeptToLca(cli::pept2lca::PeptToLca),
    #[structopt(name = "printindex")] PrintIndex(cli::printindex::PrintIndex),
    #[cfg(target_family = "unix")] #[structopt(name = "prot2kmer2lca")] ProtToKmerToLca(cli::prot2kmer2lca::ProtToKmerToLca),
    #[structopt(name = "prot2kmer")] ProtToKmer(cli::prot2kmer::ProtToKmer),
    #[structopt(name = "prot2tryp2lca")] ProtToTrypToLca(cli::prot2tryp2lca::ProtToTrypToLca),
    #[structopt(name = "prot2tryp")] ProtToTryp(cli::prot2tryp::ProtToTryp),
    #[structopt(name = "seedextend")] SeedExtend(cli::seedextend::SeedExtend),
    #[structopt(name = "snaptaxon")] SnapTaxon(cli::snaptaxon::SnapTaxon),
    #[structopt(name = "splitkmers")] SplitKmers(cli::splitkmers::SplitKmers),
    #[structopt(name = "taxa2agg")] TaxaToAgg(cli::taxa2agg::TaxaToAgg),
    #[structopt(name = "taxa2freq")] TaxaToFreq(cli::taxa2freq::TaxaToFreq),
    #[structopt(name = "taxa2tree")] TaxaToTree(cli::taxa2tree::TaxaToTree),
    #[structopt(name = "taxonomy")] Taxonomy(cli::taxonomy::Taxonomy),
    #[structopt(name = "translate")] Translate(cli::translate::Translate),
    #[structopt(name = "uniq")] Uniq(cli::uniq::Uniq),
}
