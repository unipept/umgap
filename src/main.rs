use std::collections::HashMap;
use std::io;

#[macro_use(quick_main)]
extern crate error_chain;

#[macro_use(json)]
extern crate serde_json;

use structopt::StructOpt;

use umgap::args;
use umgap::commands;
use umgap::errors::Result;
use umgap::io::fasta;
use umgap::taxon::TaxonId;

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
        args::Opt::CountRecords => countrecords(),
        args::Opt::Visualize(args) => visualize(args),
    }
});

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

fn visualize(args: args::Visualize) -> Result<()> {
    let mut taxa = HashMap::new();
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        let taxon = record.sequence[0].parse::<TaxonId>()?;
        *taxa.entry(taxon).or_insert(0) += 1;
    }

    let json = json!({
        "counts": taxa,
        "link": args.url.to_string(),
    });

    let client = reqwest::blocking::Client::new();
    let res = client
        .post("http://api.unipept.ugent.be/api/v1/taxa2tree")
        .json(&json)
        .send()
        .map_err(|err| err.to_string())?;

    print!("{}", res.text().map_err(|err| err.to_string())?);
    Ok(())
}
