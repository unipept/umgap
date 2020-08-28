//! The `umgap visualize` command.

use std::collections::HashMap;
use std::io;

use serde_json::json;

use crate::errors;
use crate::io::fasta;
use crate::taxon::TaxonId;

#[structopt(verbatim_doc_comment)]
/// Visualizes a stream of taxon IDs using the Unipept API
///
/// The `umgap visualize` command, similar to the `unipept taxa2tree` command, takes one or more
/// taxon IDs as input and returns the taxonomic tree of there taxa as output. It uses the Unipept
/// API server.
///
/// The input is given on *standard input* in a FASTA format. Per FASTA header, there should be
/// one taxon IDs. A HTML file is written to *standard output* containing a visualization of the
/// taxonomic tree containing all these taxa. With `-u`, instead a URL to the visualization hosted
/// online is printed.
///
/// ```sh
/// $ cat input.txt
/// >header1
/// 817
/// 329854
/// 1099853
/// $ umgap visualize < input.txt > output.html
/// $ umgap visualize --url < input.txt
/// https://bl.ocks.org/a686a37e1dcd43dd4ec7d467487bd6a1
/// ```
#[derive(Debug, StructOpt)]
pub struct Visualize {
    /// Host the result online and return the URL
    #[structopt(short = "u", long = "url")]
    pub url: bool,
}

/// Implements the visualize command.
pub fn visualize(args: Visualize) -> errors::Result<()> {
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

    let res = attohttpc::post("http://api.unipept.ugent.be/api/v1/taxa2tree")
        .json(&json)
        .map_err(|err| err.to_string())?
        .send()
        .map_err(|err| err.to_string())?;

    if args.url {
        let jsonres: serde_json::Value = res.json().map_err(|err| err.to_string())?;
        let gist = jsonres
            .get("gist")
            .expect("Incompatible server API")
            .as_str()
            .expect("Incompatible server API");
        println!(
            "{}",
            gist.replace("https://gist.github.com/", "https://bl.ocks.org/")
        )
    } else {
        print!("{}", res.text().map_err(|err| err.to_string())?);
    }

    Ok(())
}
