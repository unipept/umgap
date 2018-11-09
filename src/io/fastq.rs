//! Allows operations over the [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).

use std::fmt;
use std::io;
use std::io::BufRead;
use std::io::Read;
use std::iter::Peekable;

use errors;
use errors::Result;

/// Reads a FASTQ-formatted source (e.g. a file).
pub struct Reader<R: Read> {
	lines: Peekable<io::Lines<io::BufReader<R>>>,
}

impl<R: Read> Reader<R> {
	/// Creates a Reader from the given Read (e.g. a file)
	pub fn new(readable: R) -> Self {
		Reader { lines: io::BufReader::new(readable).lines().peekable() }
	}

	/// Reads the next record from the FASTQ file.
	pub fn read_record(&mut self) -> Result<Option<Record>> {
		// reading the header
		let mut header = match self.lines.next() {
			None => return Ok(None),
			Some(header) => header?,
		};
		if !header.starts_with('@') {
			bail!(errors::ErrorKind::Io(io::Error::new(
				io::ErrorKind::Other,
				"Expected @ at beginning of fastq header."
			)));
		}
		let _ = header.remove(0);

		// reading the sequence
		let mut lines = 0;
		let mut sequence = String::new();
		while self.lines
		          .peek()
		          .and_then(|line| line.as_ref().ok())
		          .map(|line| !line.starts_with('+'))
		          .unwrap_or(false)
		{
			sequence.push_str(&self.lines.next().unwrap()?);
			lines += 1;
		}

		// skipping the separator
		if self.lines
		       .next()
		       .and_then(|line| line.ok())
		       .map(|line| !line.starts_with('+'))
		       .unwrap_or(false)
		{
			bail!(errors::ErrorKind::Io(io::Error::new(
				io::ErrorKind::Other,
				"Expected a + as separator."
			)));
		}

		// reading the quality
		let mut quality = String::with_capacity(sequence.len());
		for _ in 0..lines {
			if let Some(line) = self.lines.next() {
				quality.push_str(&line?)
			} else {
				bail!(errors::ErrorKind::Io(io::Error::new(
					io::ErrorKind::Other,
					"Expected as many quality lines as sequence lines."
				)));
			}
		}

		Ok(Some(Record { header: header,
		                 sequence: sequence,
		                 quality: quality }))
	}

	/// Returns a Records struct with itself as its reader.
	pub fn records(self) -> Records<R> {
		Records { reader: self }
	}
}

/// A record as defined by the FASTQ format.
#[derive(Debug)]
pub struct Record {
	/// The FASTQ header (without the preceding '@')
	pub header: String,

	/// The actual sequence of nucleotides
	pub sequence: String,

	/// The line discribing the quality of the reads.
	/// Each character is an integer representing the estimated probability of the base being
	/// incorrect.
	pub quality: String,
}

impl fmt::Display for Record {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(
		       f,
		       "Record ({},{},{})",
		       self.header, self.sequence, self.quality
		)
	}
}

/// Convenience struct which allows for iteration (e.g. using for..in).
pub struct Records<R: Read> {
	reader: Reader<R>,
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
