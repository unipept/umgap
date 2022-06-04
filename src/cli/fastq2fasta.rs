//! The `umgap fastq2fasta` command.

use std::fs;
use std::io;
use std::path::PathBuf;

use crate::errors;
use crate::io::fasta;
use crate::io::fastq;
use crate::utils;

#[rustfmt::skip]
#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Interleaves FASTQ files into a FASTA stream
///
/// The `umgap fastq2fasta` command takes one or more FASTQ files and interleaves them into a single
/// FASTA file.
///
/// The FASTQ input files are given as command line arguments. In order, a single record is taken
/// from each of these files, and the record header and sequence are written to *standard output* in
/// FASTA format, dropping the quality scores, until any of the files runs out or records.
///
/// This command is generally used to combine two paired-end FASTQ files into a single FASTA file.
///
/// ```sh
/// $ cat input_1.fq
/// @record1/1
/// GATAGCCGTCGAGCGTCCCGACATCCATATACAGCTCGCCGCAGTGCGTCGCCCGCACGGCATCGCGCGCCGCGATGTACGCATTCAACCAGCGTCCGAG
/// +
/// AADGGGG<GI@IIKJKKKKH4EIJCHJ9:IJHKIKDIKDKGDJD@C@<>KD=;FEA:DA=I$EEED$>C@1EDE?D:CEAC;CDE:E$D$=D$EAD?AEE
/// @record2/1
/// CCCAGGTCCCCGGCATCGTCGCGGCCTCGCCCATGATCCAGCTCCACGACCAGATCCCCGTTCCCGGCGGTAAAGAGCGCGGCGTGCTCATCCTCGGAGT
/// +
/// ADDEEG@GIGIIHKCJKH@HHGKHHKHKKJJBA.GFIGK(IHKKEKECBEEEEDIKC@H<EDBJDEA36;6EE$E:G6C=E$E@CE?EE9FEE?E:F$?$
/// $ cat input_2.fq
/// @record1/2
/// CATTGTTCGCTACTTTGCGGAGCGCAATTATGCCGCGGAGATCTTCTACGTGGTGCAGCAGAAGCTGGCGGGCTGGTGCGATGCGTTGTATCGGCCCGAG
/// +
/// DDDEGGG?HIHIIKHK?@2KBHGDCJKI(JEJJKKHKKHBHKKFICEICECCFFEICCCC$E6ED$?CEEDDED$DEDCFFECEEEEFB$CCEC$6C=CA
/// @record2/2
/// GGACACGCTCTCAGGACGATGGCGCGATTGCAGGACTTGCTGGATCTCCTCCGTCGCCAAGGGGACGCGCTCGGAGTGGCTCATGGAGCAGACGAGTTCT
/// +
/// AADGGGEGIIIHIJKGCK<KD:KKHI?HIHHJKFJEKKJIGE$CKHE$EE$FEEEI=EAE8EAIKFBEE$EADEEDB$DEEDE=?B6C$C$6$A$$=BEE
/// $ umgap fastq2fasta input_1.fq input_2.fq
/// >record1/1
/// GATAGCCGTCGAGCGTCCCGACATCCATATACAGCTCGCCGCAGTGCGTCGCCCGCACGGCATCGCGCGCCGCGATGTACGCATTCAACCAGCGTCCGAG
/// >record1/2
/// CATTGTTCGCTACTTTGCGGAGCGCAATTATGCCGCGGAGATCTTCTACGTGGTGCAGCAGAAGCTGGCGGGCTGGTGCGATGCGTTGTATCGGCCCGAG
/// >record2/1
/// CCCAGGTCCCCGGCATCGTCGCGGCCTCGCCCATGATCCAGCTCCACGACCAGATCCCCGTTCCCGGCGGTAAAGAGCGCGGCGTGCTCATCCTCGGAGT
/// >record2/2
/// GGACACGCTCTCAGGACGATGGCGCGATTGCAGGACTTGCTGGATCTCCTCCGTCGCCAAGGGGACGCGCTCGGAGTGGCTCATGGAGCAGACGAGTTCT
/// ```
pub struct FastqToFasta {
    /// The input files
    #[structopt(parse(from_os_str))]
    pub input: Vec<PathBuf>,
}

/// Implements the fastq2fasta command.
pub fn fastq2fasta(args: FastqToFasta) -> errors::Result<()> {
    let handles = args
        .input
        .iter()
        .map(fs::File::open)
        .collect::<io::Result<Vec<fs::File>>>()?;
    let readers = handles
        .iter()
        .map(fastq::Reader::new)
        .map(fastq::Reader::records)
        .collect();
    let mut writer = fasta::Writer::new(io::stdout(), "", false);
    for recordzip in utils::Zip::new(readers) {
        for record in recordzip {
            let record = record?;
            writer.write_record(fasta::Record {
                header: record.header,
                sequence: vec![record.sequence],
            })?;
        }
    }
    Ok(())
}
