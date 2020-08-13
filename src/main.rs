use error_chain::quick_main;

use structopt::StructOpt;

use umgap::args;
use umgap::commands;
use umgap::errors::Result;

quick_main!(|| -> Result<()> {
    match args::Opt::from_args() {
        args::Opt::Translate(args) => commands::translate::translate(args),
        args::Opt::PeptToLca(args) => commands::pept2lca::pept2lca(args),
        args::Opt::ProtToKmerToLca(args) => commands::prot2kmer2lca::prot2kmer2lca(args),
        args::Opt::ProtToTrypToLca(args) => commands::prot2tryp2lca::prot2tryp2lca(args),
        args::Opt::TaxaToAgg(args) => commands::taxa2agg::taxa2agg(args),
        args::Opt::ProtToPept(args) => commands::prot2pept::prot2pept(args),
        args::Opt::ProtToKmer(args) => commands::prot2kmer::prot2kmer(args),
        args::Opt::Filter(args) => commands::filter::filter(args),
        args::Opt::Uniq(args) => commands::uniq::uniq(args),
        args::Opt::FastqToFasta(args) => commands::fastq2fasta::fastq2fasta(args),
        args::Opt::Taxonomy(args) => commands::taxonomy::taxonomy(args),
        args::Opt::SnapTaxon(args) => commands::snaptaxon::snaptaxon(args),
        args::Opt::SeedExtend(args) => commands::seedextend::seedextend(args),
        args::Opt::Report(args) => commands::report::report(args),
        args::Opt::BestOf(args) => commands::bestof::bestof(args),
        args::Opt::PrintIndex(args) => commands::printindex::printindex(args),
        args::Opt::SplitKmers(args) => commands::splitkmers::splitkmers(args),
        args::Opt::JoinKmers(args) => commands::joinkmers::joinkmers(args),
        args::Opt::BuildIndex(args) => commands::buildindex::buildindex(args),
        args::Opt::Visualize(args) => commands::visualize::visualize(args),
    }
});
