//! The [fastq2fasta](crate::cli::fastq2fasta::FastqToFasta) pipe.

use crate::errors;
use crate::io::fasta;
use crate::io::fastq;

/// An iterator adapter which yields a [FASTA record](crate::io::fasta::Record) for each
/// [FASTQ record](crate::io::fastq::Record) in the inputs, interleaving records from the inputs
/// if there are multiple.
pub fn pipe<I: Iterator<Item = errors::Result<fastq::Record>>>(inputs: Vec<I>) -> FastqToFasta<I> {
    FastqToFasta {
        next_input: 0,
        inputs,
    }
}

/// An iterator over FASTA records.
pub struct FastqToFasta<I: Iterator<Item = errors::Result<fastq::Record>>> {
    next_input: usize,
    inputs: Vec<I>,
}

impl<I: Iterator<Item = errors::Result<fastq::Record>>> Iterator for FastqToFasta<I> {
    type Item = errors::Result<fasta::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(fq_result) = self.inputs[self.next_input].next() {
            self.next_input = (self.next_input + 1) % self.inputs.len();
            Some(fq_result.map(|fq| fasta::Record {
                header: fq.header,
                sequence: vec![fq.sequence],
            }))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::*;

    #[test]
    fn test_empty_stays_empty() {
        let fq1 = vec![].into_iter();
        let fq2 = vec![].into_iter();
        let mut fa = pipe(vec![fq1, fq2]);
        assert_matches!(fa.next(), None);
    }

    #[test]
    fn test_stops_on_empty() {
        let fq1 = vec![Ok(fastq::Record::new("record1/1", "GATATA", ""))].into_iter();
        let fq2 = vec![].into_iter();
        let mut fa = pipe(vec![fq1, fq2]);
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record1/1", "GATATA")
        );
        assert_matches!(fa.next(), None);
    }

    #[test]
    fn test_errors_propagate() {
        let fq1 = vec![Ok(fastq::Record::new("record1/1", "GATATA", ""))].into_iter();
        let fq2 = vec![Err(ErrorKind::Test("oopsie!").into())].into_iter();
        let mut fa = pipe(vec![fq1, fq2]);
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record1/1", "GATATA")
        );
        assert_matches!(fa.next(), Some(Err(Error(ErrorKind::Test("oopsie!"), _))));
    }

    #[test]
    fn test_success() {
        let fq1 = vec![
            Ok(fastq::Record::new("record1/1", "GATATA", "")),
            Ok(fastq::Record::new("record2/1", "GATAGA", "")),
        ]
        .into_iter();
        let fq2 = vec![
            Ok(fastq::Record::new("record1/2", "TAGATA", "")),
            Ok(fastq::Record::new("record2/2", "TATAGA", "")),
        ]
        .into_iter();
        let mut fa = pipe(vec![fq1, fq2]);
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record1/1", "GATATA")
        );
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record1/2", "TAGATA")
        );
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record2/1", "GATAGA")
        );
        assert_eq!(
            fa.next().unwrap().unwrap(),
            fasta::Record::new("record2/2", "TATAGA")
        );
        assert_matches!(fa.next(), None);
    }
}
