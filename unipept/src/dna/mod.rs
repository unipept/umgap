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

impl Nucleotide {
    /// Complement of the given nucleotide.
    pub fn complement(&self) -> Self {
        match self {
            &Nucleotide::A => Nucleotide::T,
            &Nucleotide::C => Nucleotide::G,
            &Nucleotide::G => Nucleotide::C,
            &Nucleotide::T => Nucleotide::A,
            &Nucleotide::N => Nucleotide::N,
        }
    }
}

impl From<u8> for Nucleotide {
    fn from(ch: u8) -> Self {
        match ch {
            b'A' => Nucleotide::A,
            b'C' => Nucleotide::C,
            b'G' => Nucleotide::G,
            b'T' => Nucleotide::T,
            _    => Nucleotide::N,
        }
    }
}

impl Into<u8> for Nucleotide {
    fn into(self) -> u8 {
        match self {
            Nucleotide::A => b'A',
            Nucleotide::C => b'C',
            Nucleotide::G => b'G',
            Nucleotide::T => b'T',
            Nucleotide::N => b'N',
        }
    }
}

impl<'a> From<&'a u8> for Nucleotide {
    fn from(ch: &u8) -> Self {
        Nucleotide::from(*ch)
    }
}

/// A reading frame of a DNA strand.
#[derive(Clone, Copy)]
pub struct Frame<'a>(&'a [Nucleotide]);

/// A DNA strand, aka a sequence of nucleotides.
pub struct Strand(Vec<Nucleotide>);

impl<'a> From<&'a [u8]> for Strand {
    fn from(read: &[u8]) -> Self {
        Strand(read.iter().map(Nucleotide::from).collect())
    }
}

impl Strand {
    /// A reading frame of this strand, 1-indexed.
    pub fn frame<'a>(&'a self, start: usize) -> Frame<'a> {
        Frame(&self.0[start - 1..])
    }

    /// The reverse strand of this one.
    pub fn reversed(&self) -> Strand {
        Strand(self.0.iter().rev().map(Nucleotide::complement).collect())
    }
}
