//! The `umgap taxa2tree` command.

use std::collections::HashMap;
use std::io;

use serde_json::json;

use crate::errors;
use crate::io::fasta;
use crate::taxon::TaxonId;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Visualizes a stream of taxon IDs using the Unipept API
///
/// The `umgap taxa2tree` command, similar to the `unipept taxa2tree` command, takes one or more
/// taxon IDs and returns a visualized taxonomic tree of these taxa using the Unipept API server.
///
/// The input is given in a FASTA format on *standard input*. Each FASTA record contains one taxon
/// ID. A HTML file is written to *standard output* containing a visualization of the taxonomic
/// tree containing all these taxa. With `-u`, a URL to the visualization hosted online is printed
/// instead.
///
/// ```sh
/// $ cat input.txt
/// >header1
/// 817
/// 329854
/// 1099853
/// $ umgap taxa2tree < input.txt > output.html
/// $ umgap taxa2tree --url < input.txt
/// https://bl.ocks.org/a686a37e1dcd43dd4ec7d467487bd6a1
/// ```
pub struct TaxaToTree {
    /// Host the result online and return the URL
    #[structopt(short = "u", long = "url")]
    pub url: bool,
}

/// Implements the taxa2tree command.
pub fn taxa2tree(args: TaxaToTree) -> errors::Result<()> {
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
