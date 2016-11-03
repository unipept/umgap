use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Lines;
use std::io::Read;
use std::io::Write;
use std::iter::Peekable;

use errors;
use errors::Result;


const FASTA_WIDTH: usize = 70;

/// Reads a source (e.g. a file) in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
pub struct Reader<R: Read> {
    lines: Peekable<Lines<BufReader<R>>>,
    keep_lines: bool,
}


impl<R: Read> Reader<R> {
    /// Creates a Reader from the given Read (e.g. a file)
    /// keep_lines will perserve newlines in the sequence.
    pub fn new(reader: R, keep_lines: bool) -> Self {
        Reader {
            lines: BufReader::new(reader).lines().peekable(),
            keep_lines: keep_lines,
        }
    }

    /// Reads the next record from the FASTA file.
    pub fn read_record(&mut self) -> Result<Option<Record>> {
        let mut header = match self.lines.next() {
            None         => return Ok(None),
            Some(header) => try!(header),
        };

        if !header.starts_with('>') {
            return Err(errors::Error::Io(
                io::Error::new(io::ErrorKind::Other,
                               "Expected > at beginning of fasta header.")));
        }
        let _ = header.remove(0);

        let mut sequence = String::new();
        while self.lines.peek()
                        .and_then(|line| line.as_ref().ok())
                        .map(|line| !line.starts_with('>'))
                        .unwrap_or(false) {
            sequence.push_str(&try!(self.lines.next().unwrap()));
            if self.keep_lines { sequence.push('\n'); }
        }

        Ok(Some(Record { header: header, sequence: sequence.trim_right().to_string() }))
    }

    /// Returns a Records struct with itself as its reader.
    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

/// A record as defined by the FASTA format.
#[derive(Debug)]
pub struct Record {
    /// The record header (without the preceding '<')
    pub header: String,

    /// The actual sequence of nucleotides
    pub sequence: String,
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
pub struct Writer<W: Write> {
    buffer: io::BufWriter<W>,
    wrap: bool,
}

impl<W: Write> Writer<W> {
    pub fn new(writer: W, wrap: bool) -> Self {
        Writer {
            buffer: io::BufWriter::new(writer),
            wrap: wrap,
        }
    }

    /// Convenience method, see [write_record_ref](#method.write_record_ref).
    pub fn write_record(&mut self, record: Record) -> Result<()> {
        self.write_record_ref(&record)
    }

    /// Writes a Record to the Write, in FASTA format
    pub fn write_record_ref(&mut self, record: &Record) -> Result<()> {
        try!(write!(self.buffer, ">{}\n", record.header));
        if !self.wrap {
            try!(self.buffer.write(record.sequence.as_bytes()));
            try!(self.buffer.write_all(&[b'\n']));
        } else {
            for subseq in record.sequence.as_bytes().chunks(FASTA_WIDTH) {
                try!(self.buffer.write_all(subseq));
                try!(self.buffer.write_all(&[b'\n']));
            }
        }
        Ok(())
    }
}

