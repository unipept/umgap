use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Lines;
use std::io::Read;
use std::io::Write;

use errors;
use errors::Result;


const FASTA_WIDTH: usize = 70;

pub struct Reader<R: Read> {
    lines: Lines<BufReader<R>>,
}


impl<R: Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            lines: BufReader::new(reader).lines(),
        }
    }

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


pub struct Record {
    pub header: String,
    pub sequence: String,
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


pub struct Writer<W: Write> {
    buffer: io::BufWriter<W>,
}

impl<W: Write> Writer<W> {
    pub fn new(writer: W) -> Self {
        Writer {
            buffer: io::BufWriter::new(writer),
        }
    }

    pub fn write_record(&mut self, record: Record) -> Result<()> {
        self.write_record_ref(&record)
    }

    pub fn write_record_ref(&mut self, record: &Record) -> Result<()> {
        try!(write!(self.buffer, ">{}\n", record.header));
        for subseq in record.sequence.as_bytes().chunks(FASTA_WIDTH) {
            try!(self.buffer.write_all(subseq));
            try!(self.buffer.write_all(&[b'\n']));
        }
        Ok(())
    }
}

