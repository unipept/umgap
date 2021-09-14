//! The `umgap translate` command

use std::fmt;
use std::io;
use std::str::FromStr;

use crate::dna::translation::TranslationTable;
use crate::dna::Strand;
use crate::errors;
use crate::io::fasta;

#[derive(Debug, StructOpt)]
#[structopt(verbatim_doc_comment)]
/// Translates a FASTA stream of DNA reads to peptides
///
/// The `umgap translate` command takes one or more DNA sequences and translates them into amino
/// acid sequences.
///
/// The DNA sequences are expected in a FASTA format on *standard input*.
///
/// ```sh
/// $ cat input.fa
/// >header1
/// GATTACAAA
/// $ umgap translate -f1 < input.fa
/// >header1
/// DYK
/// ```
///
/// The `-f` flag allows you to add reading frames to the translation. If you want to translate
/// multiple frames and care to keep them apart, the `-n` flag adds the name of the frame to the end
/// of the FASTA header:
///
/// ```sh
/// $ umgap translate -f1 -f1R < input.fa
/// >header1|1
/// DYK
/// >header1|1R
/// FVI
/// ```
///
/// The `-a` flag can be used as a shortcut to translate all 6 reading frames.
///
/// With the `-t` flag, you can select a specific translation table, for instance `-t11` for the
/// bacterial, archaeal and plant plastid code.
pub struct Translate {
    /// Replace each start-codon with methionine
    #[structopt(short = "m", long = "methionine")]
    pub methionine: bool,

    /// Read and output all six frames
    #[structopt(short = "a", long = "all-frames", conflicts_with = "frame")]
    pub all_frames: bool,

    /// Adds a reading frame (1, 2, 3, 1R, 2R or 3R)
    #[structopt(
        short = "f",
        long = "frame",
        possible_values = &Frame::variants()
    )]
    pub frames: Vec<Frame>,

    /// Append a bar (|) and the name of the frame to the fasta header
    #[structopt(short = "n", long = "append-name")]
    pub append_name: bool,

    /// Translation table to use
    #[structopt(short = "t", long = "table", default_value = "1")]
    pub table: String,

    /// Instead of normal use, print the selected table and exit
    #[structopt(short = "s", long = "show-table")]
    pub show_table: bool,
}

/// Implements the translate command
pub fn translate(args: Translate) -> errors::Result<()> {
    // Parsing the table
    let table = args.table.parse::<&TranslationTable>()?;

    // Which frames to do
    let frames = if args.all_frames {
        vec![
            Frame::Forward1,
            Frame::Forward2,
            Frame::Forward3,
            Frame::Reverse1,
            Frame::Reverse2,
            Frame::Reverse3,
        ]
    } else {
        args.frames
    };

    // Split on show_tables
    if args.show_table {
        table.print();
    } else {
        let mut writer = fasta::Writer::new(io::stdout(), "", false);

        // Parsing the frames
        let frames = frames
            .iter()
            .map(|&frame| match frame {
                Frame::Forward1 => (frame, 1, false),
                Frame::Forward2 => (frame, 2, false),
                Frame::Forward3 => (frame, 3, false),
                Frame::Reverse1 => (frame, 1, true),
                Frame::Reverse2 => (frame, 2, true),
                Frame::Reverse3 => (frame, 3, true),
            })
            .collect::<Vec<(Frame, usize, bool)>>();

        for record in fasta::Reader::new(io::stdin(), true).records() {
            let fasta::Record { header, sequence } = record?;

            let forward = Strand::from(&sequence);
            let reverse = forward.reversed();
            for &(name, frame, reversed) in &frames {
                let strand = if reversed { &reverse } else { &forward };
                writer.write_record(fasta::Record {
                    header: if !args.append_name {
                        header.clone()
                    } else {
                        header.clone() + "|" + &name.to_string()
                    },
                    sequence: vec![String::from_utf8(
                        table.translate_frame(args.methionine, strand.frame(frame)),
                    )
                    .unwrap()],
                })?;
            }
        }
    }
    Ok(())
}

/// A reading frame
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy)]
pub enum Frame {
    Forward1,
    Forward2,
    Forward3,
    Reverse1,
    Reverse2,
    Reverse3,
}

static FRAMES: &[&str] = &["1", "2", "3", "1R", "2R", "3R"];
impl Frame {
    fn variants() -> &'static [&'static str] {
        FRAMES
    }
}

impl FromStr for Frame {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "1" => Ok(Frame::Forward1),
            "2" => Ok(Frame::Forward2),
            "3" => Ok(Frame::Forward3),
            "1R" => Ok(Frame::Reverse1),
            "2R" => Ok(Frame::Reverse2),
            "3R" => Ok(Frame::Reverse3),
            _ => Err(ErrorKind::ParseFrameError(s.to_string()).into()),
        }
    }
}

impl fmt::Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match *self {
                Frame::Forward1 => "1",
                Frame::Forward2 => "2",
                Frame::Forward3 => "3",
                Frame::Reverse1 => "1R",
                Frame::Reverse2 => "2R",
                Frame::Reverse3 => "3R",
            }
        )
    }
}

error_chain! {
    errors {
        /// Unparseable Frame
        ParseFrameError(frame: String) {
            description("Unparseable frame")
            display("Unparseable frame: {}", frame)
        }
    }
}
