use error_chain::quick_main;

use structopt::StructOpt;

use umgap::commands;
use umgap::errors::Result;

quick_main!(|| -> Result<()> {
    match Opt::from_args() {
        Opt::Translate(args) => commands::translate::translate(args),
        Opt::PeptToLca(args) => commands::pept2lca::pept2lca(args),
        Opt::ProtToKmerToLca(args) => commands::prot2kmer2lca::prot2kmer2lca(args),
        Opt::ProtToTrypToLca(args) => commands::prot2tryp2lca::prot2tryp2lca(args),
        Opt::TaxaToAgg(args) => commands::taxa2agg::taxa2agg(args),
        Opt::ProtToPept(args) => commands::prot2pept::prot2pept(args),
        Opt::ProtToKmer(args) => commands::prot2kmer::prot2kmer(args),
        Opt::Filter(args) => commands::filter::filter(args),
        Opt::Uniq(args) => commands::uniq::uniq(args),
        Opt::FastqToFasta(args) => commands::fastq2fasta::fastq2fasta(args),
        Opt::Taxonomy(args) => commands::taxonomy::taxonomy(args),
        Opt::SnapTaxon(args) => commands::snaptaxon::snaptaxon(args),
        Opt::SeedExtend(args) => commands::seedextend::seedextend(args),
        Opt::Report(args) => commands::report::report(args),
        Opt::BestOf(args) => commands::bestof::bestof(args),
        Opt::PrintIndex(args) => commands::printindex::printindex(args),
        Opt::SplitKmers(args) => commands::splitkmers::splitkmers(args),
        Opt::JoinKmers(args) => commands::joinkmers::joinkmers(args),
        Opt::BuildIndex(args) => commands::buildindex::buildindex(args),
        Opt::Visualize(args) => commands::visualize::visualize(args),
    }
});

/// The Options enum for UMGAP arguments
#[derive(Debug, StructOpt)]
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
    ProtToTrypToLca(commands::prot2tryp2lca::ProtToTrypToLca),

    /// Aggregates taxa to a single taxon.
    #[structopt(name = "taxa2agg")]
    TaxaToAgg(commands::taxa2agg::TaxaToAgg),

    /// Splits each protein sequence in a FASTA format into a list of (tryptic) peptides.
    #[structopt(name = "prot2pept")]
    ProtToPept(commands::prot2pept::ProtToPept),

    /// Pick the frame with the most none-root hits.
    #[structopt(name = "bestof")]
    BestOf(commands::bestof::BestOf),

    /// Count and report on a list of taxon ids.
    #[structopt(name = "report")]
    Report(commands::report::Report),

    /// Seed and extend.
    #[structopt(name = "seedextend")]
    SeedExtend(commands::seedextend::SeedExtend),

    /// Show taxonomy info
    #[structopt(name = "taxonomy")]
    Taxonomy(commands::taxonomy::Taxonomy),

    /// Snap taxa to a specified rank or one of the specified taxa.
    #[structopt(name = "snaptaxon")]
    SnapTaxon(commands::snaptaxon::SnapTaxon),

    /// Interleaves a number of FASTQ files into a single FASTA output.
    #[structopt(name = "fastq2fasta")]
    FastqToFasta(commands::fastq2fasta::FastqToFasta),

    /// Filter peptides in a FASTA format based on specific criteria.
    #[structopt(name = "filter")]
    Filter(commands::filter::Filter),

    /// Concatenates the data strings of all consecutive FASTA entries with the same header.
    #[structopt(name = "uniq")]
    Uniq(commands::uniq::Uniq),

    /// Splits each protein sequence in a FASTA format into a list of kmers.
    #[structopt(name = "prot2kmer")]
    ProtToKmer(commands::prot2kmer::ProtToKmer),

    /// Print the values in an FST index to stdout.
    #[structopt(name = "printindex")]
    PrintIndex(commands::printindex::PrintIndex),

    /// Splits each protein sequence in a CSV format into a list of kmers.
    #[structopt(name = "splitkmers")]
    SplitKmers(commands::splitkmers::SplitKmers),

    /// Groups a CSV by equal first fields (Kmers) and aggregates the second fields (taxon ids).
    #[structopt(name = "joinkmers")]
    JoinKmers(commands::joinkmers::JoinKmers),

    /// Write an FST index of stdin on stdout.
    #[structopt(name = "buildindex")]
    BuildIndex(commands::buildindex::BuildIndex),

    /// Visualizes the given list of taxons using the Unipept API
    #[structopt(name = "visualize")]
    Visualize(commands::visualize::Visualize),
}
