use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use errors;
use errors::Result;


pub struct Reader<R: io::Read> {
    reader: BufReader<R>,
}


impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            reader: BufReader::new(reader),
        }
    }

    pub fn read_record(&mut self) -> Result<Option<Record>> {
        let mut header = String::new();
        let bytes_read = try!(self.reader.read_line(&mut header));

        if bytes_read == 0 {
            return Ok(None);
        }

        if !header.starts_with('>') {
            return Err(errors::Error::Io(
                io::Error::new(io::ErrorKind::Other,
                               "Expected > at beginning of fasta header.")));
        }

        let mut sequence = String::new();
        let bytes_read = try!(self.reader.read_line(&mut sequence));

        if bytes_read == 0 || sequence.starts_with('>') {
            return Err(errors::Error::Io(
                io::Error::new(io::ErrorKind::Other,
                               "Encountered empty sequence.")));
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
