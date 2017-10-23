//! Allows operations over the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).


use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Lines;
use std::io::Read;
use std::io::Write;
use std::iter::Peekable;

use regex::Regex;

use errors;
use errors::Result;


const FASTA_WIDTH: usize = 70;

/// Reads a FASTA-formatted source (e.g. a file).
pub struct Reader<R: Read> {
    lines: Peekable<Lines<BufReader<R>>>,
    separator: Option<Regex>,
    wrapped: bool,
}


impl<R: Read> Reader<R> {
    /// Creates a Reader from the given Read (e.g. a file). Passing a separator will split the
    /// fasta entry on this separator.
    pub fn new(reader: R, separator: Option<Regex>, wrapped: bool) -> Self {
        Reader {
            lines: BufReader::new(reader).lines().peekable(),
            separator: separator,
            wrapped: wrapped,
        }
    }

    /// Reads the next record from the FASTA file.
    pub fn read_record(&mut self) -> Result<Option<Record>> {
        let mut header = match self.lines.next() {
            None         => return Ok(None),
            Some(header) => header?,
        };

        if !header.starts_with('>') {
            bail!(errors::ErrorKind::Io(
                io::Error::new(io::ErrorKind::Other,
                               "Expected > at beginning of fasta header.")));
        }
        let _ = header.remove(0);

        let mut sequence = String::new();
        while self.lines.peek()
                        .and_then(|line| line.as_ref().ok())
                        .map(|line| !line.starts_with('>'))
                        .unwrap_or(false) {
            sequence.push_str(&self.lines.next().unwrap()?);
            if !self.wrapped { sequence.push('\n') }
        }

        Ok(Some(Record {
            header: header,
            sequence: self.separator.as_ref()
                                    .map(|re| re.split(sequence.trim_right())
                                                .map(str::to_string)
                                                .collect())
                                    .unwrap_or(vec![sequence]),
        }))
    }

    /// Returns a Records struct with itself as its reader.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

/// A record extending the FASTA format. A single header can be followed by multiple
/// sequences/items.
#[derive(Debug)]
pub struct Record {
    /// The record header (without the preceding '>')
    pub header: String,

    /// The actual sequence of nucleotides
    pub sequence: Vec<String>,
}

/// Convenience struct which allows for iteration (e.g. using for..in).
pub struct Records<R: Read> {
    reader: Reader<R>,
}

impl<R: Read> Iterator for Records<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.reader.read_record() {
            Ok(None)         => None,
            Ok(Some(record)) => Some(Ok(record)),
            Err(err)         => Some(Err(err)),
        }
    }
}

/// Writes to a file in the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
pub struct Writer<'a, W: Write> {
    buffer: io::BufWriter<W>,
    separator: &'a str,
    wrap: bool,
}

impl<'a, W: Write> Writer<'a, W> {
    /// Constructs a writer from the specified Write. Items are printed
    /// separated by separator.
    pub fn new(write: W, separator: &'a str, wrap: bool) -> Self {
        Writer {
            buffer: io::BufWriter::new(write),
            separator: separator,
            wrap: wrap,
        }
    }

    /// Convenience method, see [write_record_ref](#method.write_record_ref).
    pub fn write_record(&mut self, record: Record) -> Result<()> {
        self.write_record_ref(&record)
    }

    /// Writes a Record to the Write, in FASTA format
    pub fn write_record_ref(&mut self, record: &Record) -> Result<()> {
        write!(self.buffer, ">{}\n", record.header)?;
        let sequence = record.sequence.join(self.separator);
        if !self.wrap {
            self.buffer.write(sequence.as_bytes())?;
            self.buffer.write_all(&[b'\n'])?;
        } else {
            for subseq in sequence.as_bytes().chunks(FASTA_WIDTH) {
                self.buffer.write_all(subseq)?;
                self.buffer.write_all(&[b'\n'])?;
            }
        }
        Ok(())
    }
}

