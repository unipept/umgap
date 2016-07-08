
use std::io;
use std::io::Read;
use std::io::BufRead;
use std::iter::Peekable;
use std::fmt;

use errors;
use errors::Result;

pub struct Reader<R: Read> {
    lines: Peekable<io::Lines<io::BufReader<R>>>,
}

impl<R: Read> Reader<R> {
    pub fn new(readable: R) -> Self {
        Reader {
            lines: io::BufReader::new(readable).lines().peekable(),
        }
    }

    pub fn read_record(&mut self) -> Result<Option<Record>> {
        // reading the header
        let mut header = match self.lines.next() {
            None         => return Ok(None),
            Some(header) => try!(header)
        };
        if !header.starts_with('@') {
            return Err(errors::Error::Io(io::Error::new(
                io::ErrorKind::Other,
                "Expected @ at beginning of fasta header."
            )));
        }
        let _ = header.remove(0);

        // reading the sequence
        let mut lines = 0;
        let mut sequence = String::new();
        while self.lines.peek()
                  .and_then(|line| line.as_ref().ok())
                  .map(|line| !line.starts_with('+'))
                  .unwrap_or(false) {
            sequence.push_str(&try!(self.lines.next().unwrap()));
            lines += 1;
        }

        // skipping the separator
        if self.lines.next()
               .and_then(|line| line.ok())
               .map(|line| !line.starts_with('+'))
               .unwrap_or(false) {
            return Err(errors::Error::Io(io::Error::new(
                io::ErrorKind::Other,
                "Expected a + as separator."
            )));
        }

        // reading the quality
        let mut quality = String::with_capacity(sequence.len());
        for _ in 0..lines {
            if let Some(line) = self.lines.next() {
                quality.push_str(&try!(line))
            } else {
                return Err(errors::Error::Io(io::Error::new(
                    io::ErrorKind::Other,
                    "Expected as many quality lines as sequence lines."
                )));
            }
        }

        Ok(Some(Record {
            header: header,
            sequence: sequence,
            quality: quality,
        }))
    }

    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

#[derive(Debug)]
pub struct Record {
    pub header: String,
    pub sequence: String,
    pub quality: String,
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Record ({},{},{})", self.header, self.sequence, self.quality)
    }
}

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

