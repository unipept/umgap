//! Module for DNA related code.

pub mod translation;

/// A Deoxyribose nucleotide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nucleotide {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
    /// Any Nucleotide or None
    N,
}
use self::Nucleotide::*;

impl Nucleotide {
    /// Complement of the given nucleotide.
    pub fn complement(&self) -> Self {
        match self {
            &A => T,
            &C => G,
            &G => C,
            &T => A,
            &N => N,
        }
    }
}

impl From<u8> for Nucleotide {
    fn from(ch: u8) -> Self {
        match ch {
            b'A' => A,
            b'C' => C,
            b'G' => G,
            b'T' => T,
            _ => N,
        }
    }
}

impl Into<u8> for Nucleotide {
    fn into(self) -> u8 {
        match self {
            A => b'A',
            C => b'C',
            G => b'G',
            T => b'T',
            N => b'N',
        }
    }
}

impl<'a> From<&'a u8> for Nucleotide {
    fn from(ch: &u8) -> Self {
        Nucleotide::from(*ch)
    }
}

/// A reading frame of a DNA strand.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame<'a>(&'a [Nucleotide]);

/// A DNA strand, aka a sequence of nucleotides.
#[derive(Debug, PartialEq)]
pub struct Strand(Vec<Nucleotide>);

impl<'a> From<&'a [u8]> for Strand {
    fn from(read: &[u8]) -> Self {
        Strand(read.iter().map(Nucleotide::from).collect())
    }
}

impl<'a> From<&'a Vec<String>> for Strand {
    fn from(lines: &Vec<String>) -> Self {
        Strand(
            lines
                .iter()
                .flat_map(|s| s.as_bytes())
                .map(Nucleotide::from)
                .collect(),
        )
    }
}

impl Strand {
    /// A reading frame of this strand, 1-indexed.
    pub fn frame<'a>(&'a self, start: usize) -> Frame<'a> {
        Frame(if self.0.len() > start - 1 {
            &self.0[start - 1..]
        } else {
            &[]
        })
    }

    /// The reverse strand of this one.
    pub fn reversed(&self) -> Strand {
        Strand(self.0.iter().rev().map(Nucleotide::complement).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! strand {
        ( $( $x:path )* ) => (Strand(vec![$($x),*]))
    }

    #[test]
    fn test_strand_from() {
        assert_eq!(
            strand![A C G T N T C G A],
            Strand::from("ACGT*TCGA".as_bytes())
        );
    }

    #[test]
    fn test_strand_reversed() {
        assert_eq!(
            strand![A C G T N T G C A],
            strand![T G C A N A C G T].reversed()
        );
    }

    #[test]
    fn test_strand_frames() {
        assert_eq!(Frame(&[A, C, G, T]), strand![A C G T].frame(1));
        assert_eq!(Frame(&[C, G, T]), strand![A C G T].frame(2));
        assert_eq!(Frame(&[G, T]), strand![A C G T].frame(3));
    }
}
