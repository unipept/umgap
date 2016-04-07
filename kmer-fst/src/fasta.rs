use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::Lines;
use std::io::BufReader;
use std::path::Path;

use errors;
use errors::Result;


pub struct Reader<R: io::Read> {
    lines: Lines<BufReader<R>>,
}


impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            lines: BufReader::new(reader).lines(),
        }
    }

    pub fn read_record(&mut self) -> Result<Option<Record>> {
        let header = match self.lines.next() {
            None         => return Ok(None),
            Some(header) => try!(header),
        };

        if !header.starts_with('>') {
            return Err(errors::Error::Io(
                io::Error::new(io::ErrorKind::Other,
                               "Expected > at beginning of fasta header.")));
        }

        let sequence = match self.lines.next() {
            None => return Err(errors::Error::Io(io::Error::new(io::ErrorKind::Other, "Encountered empty sequence at end of file."))),
            Some(sequence) => try!(sequence),
        };

        if sequence.starts_with('>') {
            return Err(errors::Error::Io(io::Error::new(io::ErrorKind::Other, "Encountered empty sequence.")));
        }

        Ok(Some(Record { header: header, sequence: sequence }))
    }

    pub fn records(self) -> Records<R> {
        Records { reader: self }
    }
}

impl Reader<File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        File::open(path).map(Reader::new)
    }
}


pub struct Record {
    pub header: String,
    pub sequence: String,
}


pub struct Records<R: io::Read> {
    reader: Reader<R>,
}


impl<R: io::Read> Iterator for Records<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.reader.read_record() {
            Ok(None)         => None,
            Ok(Some(record)) => Some(Ok(record)),
            Err(err)         => Some(Err(err)),
        }
    }
}
