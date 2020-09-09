//! The `umgap bestof` command.

use std::io;

use crate::errors;
use crate::io::fasta;
use crate::taxon::TaxonId;

#[structopt(verbatim_doc_comment)]
/// Selects the best read of every fixed size group
///
/// The `umgap bestof` command takes groups of taxon IDs as input and outputs for each group the
/// taxon ID with the most non-root identifications.
///
/// The input is given in FASTA format on *standard input*. Per FASTA header, there should be
/// multiple numbers (taxon IDs). Per 6 FASTA records (or whichever number you specify with `-f`),
/// the best record is selected and written to *standard output*. If the input is a series of
/// identified taxon IDs for each of the 6 translations of a read, the output will most likely come
/// from the actual coding frame.
///
/// ```sh
/// $ cat dna.fa
/// >header1
/// CGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCCGGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCC
/// $ umgap translate -n -a < dna.fa | umgap prot2kmer2lca 9mer.index | tee input.fa
/// >header1|1
/// 9606 9606 2759 9606 9606 9606 9606 9606 9606 9606 8287
/// >header1|2
/// 2026807 888268 186802 1598 1883
/// >header1|3
/// 1883
/// >header1|1R
/// 27342 2759 155619 1133106 38033 2
/// >header1|2R
/// >header1|3R
/// 2951
/// $ umgap bestof < input.fa
/// >header1|1
/// 9606 9606 2759 9606 9606 9606 9606 9606 9606 9606 8287
/// ```
///
/// Taxon IDs are separated by newlines in the actual output, but are separated by spaces in this
/// example.
#[derive(Debug, StructOpt)]
pub struct BestOf {
    /// The number of frames of which to pick the best
    #[structopt(short = "f", long = "frames", default_value = "6")]
    pub frames: usize,
}

/// Implements the bestof command.
pub fn bestof(args: BestOf) -> errors::Result<()> {
    let mut writer = fasta::Writer::new(io::stdout(), "\n", false);

    // Combine frames and process them
    let mut chunk = Vec::with_capacity(args.frames);
    for record in fasta::Reader::new(io::stdin(), false).records() {
        let record = record?;
        if chunk.len() < args.frames - 1 {
            chunk.push(record);
        } else {
            // process chunk
            writer.write_record_ref(
                chunk
                    .iter()
                    .max_by_key(|&rec| {
                        rec.sequence
                            .iter()
                            .map(|tid| tid.parse::<TaxonId>().unwrap_or(0))
                            .filter(|&s| s != 0 && s != 1)
                            .count()
                    })
                    .unwrap(),
            )?;
            chunk.clear();
        }
    }
    Ok(())
}
