
use std::process;
use std::io;
use std::io::Write;
use std::collections::HashMap;
use std::str;

extern crate clap;
use clap::App;

extern crate unipept;
use unipept::{PKG_NAME, PKG_VERSION, PKG_AUTHORS};
use unipept::io::fasta;

const ABOUT: &'static str = "
Translates RNA into Amino Acid Sequences
";

fn main() {
    let _       = App::new(PKG_NAME.to_string() + " rna2aa")
                      .version(PKG_VERSION)
                      .author(PKG_AUTHORS.split(':').next().unwrap_or("unknown"))
                      .about(ABOUT)
                      .get_matches();
    main_result().unwrap_or_else(|err| {
        writeln!(&mut io::stderr(), "{}", err).unwrap();
        process::exit(1);
    });
}

fn main_result() -> Result<(), String> {
    let table = TranslationTable::new();
    let mut writer = fasta::Writer::new(io::stdout(), "", false);
    for record in fasta::Reader::new(io::stdin(), None, true).records() {
        let fasta::Record { header, sequence } = try!(record.map_err(|err| format!("Something went wrong during the reading: {}", err)));
        try!(writer.write_record(fasta::Record {
            header: header,
            sequence: sequence[0].as_bytes()
                                 .chunks(3)
                                 .map(str::from_utf8)
                                 .map(|res| res.map(|c| table.translate(c)).unwrap_or("").to_string())
                                 .collect()
        }).map_err(|err| err.to_string()));
    }

    Ok(())
}

struct TranslationTable {
    table: HashMap<&'static str, &'static str>
}

impl TranslationTable {
    fn new() -> Self {
        let mut table = HashMap::new();
        table.insert("TTT", "F");
        table.insert("TTC", "F");
        table.insert("TTA", "L");
        table.insert("TTG", "L");
        table.insert("TCT", "S");
        table.insert("TCC", "S");
        table.insert("TCA", "S");
        table.insert("TCG", "S");
        table.insert("TAT", "Y");
        table.insert("TAC", "Y");
        table.insert("TAA", "*");
        table.insert("TAG", "*");
        table.insert("TGT", "C");
        table.insert("TGC", "C");
        table.insert("TGA", "*");
        table.insert("TGG", "W");

        table.insert("CTT", "L");
        table.insert("CTC", "L");
        table.insert("CTA", "L");
        table.insert("CTG", "L");
        table.insert("CCT", "P");
        table.insert("CCC", "P");
        table.insert("CCA", "P");
        table.insert("CCG", "P");
        table.insert("CAT", "H");
        table.insert("CAC", "H");
        table.insert("CAA", "Q");
        table.insert("CAG", "Q");
        table.insert("CGT", "R");
        table.insert("CGC", "R");
        table.insert("CGA", "R");
        table.insert("CGG", "R");

        table.insert("ATT", "I");
        table.insert("ATC", "I");
        table.insert("ATA", "I");
        table.insert("ATG", "M");
        table.insert("ACT", "T");
        table.insert("ACC", "T");
        table.insert("ACA", "T");
        table.insert("ACG", "T");
        table.insert("AAT", "N");
        table.insert("AAC", "N");
        table.insert("AAA", "K");
        table.insert("AAG", "K");
        table.insert("AGT", "S");
        table.insert("AGC", "S");
        table.insert("AGA", "R");
        table.insert("AGG", "R");

        table.insert("GTT", "V");
        table.insert("GTC", "V");
        table.insert("GTA", "V");
        table.insert("GTG", "V");
        table.insert("GCT", "A");
        table.insert("GCC", "A");
        table.insert("GCA", "A");
        table.insert("GCG", "A");
        table.insert("GAT", "D");
        table.insert("GAC", "D");
        table.insert("GAA", "E");
        table.insert("GAG", "E");
        table.insert("GGT", "G");
        table.insert("GGC", "G");
        table.insert("GGA", "G");
        table.insert("GGG", "G");

        TranslationTable { table }
    }

    fn translate(&self, codon: &str) -> &str {
        self.table.get(codon).unwrap_or(&"-")
    }
}
