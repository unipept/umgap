//! Allows operations over the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Lines;
use std::io::Read;
use std::io::Write;
use std::iter::Peekable;

use crate::errors;
use crate::errors::Result;

const BUFFER_SIZE: usize = 10_000_000; // 10MB
const FASTA_WIDTH: usize = 70;

/// Reads a FASTA-formatted source (e.g. a file).
pub struct Reader<R: Read> {
    lines: Peekable<Lines<BufReader<R>>>,
    unwrap: bool,
}

impl<R: Read> Reader<R> {
    /// Creates a Reader from the given Read (e.g. a file).
    /// When unwrap is `true`, the reader will unwrap the input sequences by
    /// removing newlines between the sequence lines
    /// (each [Record.sequence](struct.Record.html)) will only have one item.
    /// When unwrap is `false`, each record line will be a new item in
    /// the [Record.sequence](struct.Record.html) vec.
    pub fn new(reader: R, unwrap: bool) -> Self {
        let lines = BufReader::with_capacity(BUFFER_SIZE, reader)
            .lines()
            .peekable();
        Reader { lines, unwrap }
    }

    /// Reads the next record from the FASTA file.
    pub fn read_record(&mut self) -> Result<Option<Record>> {
        let mut header = match self.lines.next() {
            None => return Ok(None),
            Some(header) => header?,
        };

        if !header.starts_with('>') {
            bail!(errors::ErrorKind::Io(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at beginning of fasta header."
            )));
        }
        let _ = header.remove(0);

        let mut sequence = Vec::new();
        while self
            .lines
            .peek()
            .and_then(|line| line.as_ref().ok())
            .map(|line| !line.starts_with('>'))
            .unwrap_or(false)
        {
            sequence.push(self.lines.next().unwrap()?);
        }
        if self.unwrap {
            sequence = vec![sequence.concat()];
        }

        Ok(Some(Record { header, sequence }))
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

impl<R: Read> Records<R> {
    /// Convert this to a chunked version which processes the records in chunks
    /// of the given size.
    pub fn chunked(self, size: usize) -> ChunkedRecords<R> {
        ChunkedRecords {
            reader: self.reader,
            chunk_size: size,
        }
    }
}

impl<R: Read> Iterator for Records<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.reader.read_record() {
            Ok(None) => None,
            Ok(Some(record)) => Some(Ok(record)),
            Err(err) => Some(Err(err)),
        }
    }
}

/// Allows for iteration over records in chunks
pub struct ChunkedRecords<R: Read> {
    reader: Reader<R>,
    chunk_size: usize,
}

impl<R: Read> Iterator for ChunkedRecords<R> {
    type Item = Result<Vec<Record>>;

    fn next(&mut self) -> Option<Result<Vec<Record>>> {
        let mut chunk = Vec::with_capacity(self.chunk_size);
        while chunk.len() < self.chunk_size {
            match self.reader.read_record() {
                Ok(Some(result)) => chunk.push(result),
                Ok(None) => break,
                Err(err) => return Some(Err(err)),
            }
        }
        if chunk.is_empty() {
            None
        } else {
            Some(Ok(chunk))
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
            separator,
            wrap,
        }
    }

    /// Convenience method, see [write_record_ref](#method.write_record_ref).
    pub fn write_record(&mut self, record: Record) -> Result<()> {
        self.write_record_ref(&record)
    }

    /// Writes a Record to the Write, in FASTA format
    pub fn write_record_ref(&mut self, record: &Record) -> Result<()> {
        write!(self.buffer, ">{}", record.header)?;
        let sequence = record.sequence.join(self.separator);
        if !self.wrap {
            self.buffer.write_all(&[b'\n'])?;
            self.buffer.write_all(sequence.as_bytes())?;
        } else {
            for subseq in sequence.as_bytes().chunks(FASTA_WIDTH) {
                self.buffer.write_all(&[b'\n'])?;
                self.buffer.write_all(subseq)?;
            }
        }
        if !sequence.is_empty() {
            self.buffer.write_all(&[b'\n'])?;
        }
        Ok(())
    }
}
