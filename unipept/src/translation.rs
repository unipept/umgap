//! Defines translation tables.

use std;
use std::str;
use std::collections::HashMap;

const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

// TODO use a nucleotide enum instead of u8
// #[derive(Debug)]
// enum Nucleotide { T, C, A, G }
// fn main() {
//     let s: String = format!("{:?}", Nucleotide::T);
//     println!("{}", s);
// }

/// Invert the given nucleotide.
pub fn invert(base: u8) -> u8 {
    match base {
        b'T' => b'A',
        b'C' => b'G',
        b'A' => b'T',
        b'G' => b'C',
        _    => unreachable!()
    }
}

lazy_static! {
    static ref TABLES: [Option<TranslationTable>; 23] = [
        Some(TranslationTable::new("universal", 1,
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "---M---------------M---------------M----------------------------")),
        Some(TranslationTable::new("vertebrate_mitochondrial", 2,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
            "--------------------------------MMMM---------------M------------")),
        Some(TranslationTable::new("yeast_mitochondrial", 3,
            "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "----------------------------------MM----------------------------")),
        Some(TranslationTable::new("mold_mitochondrial", 4,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "--MM---------------M------------MMMM---------------M------------")),
        Some(TranslationTable::new("invertebrate_mitochondrial", 5,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
            "---M----------------------------MMMM---------------M------------")),
        Some(TranslationTable::new("ciliate_nuclear", 6,
            "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        None, None,
        Some(TranslationTable::new("echinoderm_mitochondrial", 9,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            "-----------------------------------M---------------M------------")),
        Some(TranslationTable::new("euplotid_nuclear", 10,
            "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        Some(TranslationTable::new("bacterial", 11,
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "---M---------------M------------MMMM---------------M------------")),
        Some(TranslationTable::new("alternative_yeast_nuclear", 12,
            "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-------------------M---------------M----------------------------")),
        Some(TranslationTable::new("ascidian_mitochondrial", 13,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
            "---M------------------------------MM---------------M------------")),
        Some(TranslationTable::new("flatworm_mitochondrial", 14,
            "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        Some(TranslationTable::new("blepharisma_macronuclear", 15,
            "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        Some(TranslationTable::new("chlorophycean_mitochondrial", 16,
            "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        None, None, None, None,
        Some(TranslationTable::new("trematode_mitochondrial", 21,
            "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            "-----------------------------------M---------------M------------")),
        Some(TranslationTable::new("scenedesmus_mitochondrial", 22,
            "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "-----------------------------------M----------------------------")),
        Some(TranslationTable::new("thraustochytrium_mitochondrial", 23,
            "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            "--------------------------------M--M---------------M------------")),
    ];
}

/// A translation table.
#[derive(Debug, PartialEq)]
pub struct TranslationTable {
    name: String,
    index: u32,
    table: HashMap<Codon, (bool, u8)>
}

impl TranslationTable {
    fn new(name: &str, index: u32, aas: &str, starts: &str) -> Self {
        let name = name.to_string();
        let mut table = HashMap::new();
        for (codon, (aa, start)) in CodonIterator::new().zip(aas.bytes().zip(starts.bytes())) {
            table.insert(codon, (start == b'M', aa));
        }
        TranslationTable { name, index, table }
    }

    /// Translate the given codon to an AA. Translate each start codon to methionine if asked.
    pub fn translate(&self, methionine: bool, codon: &Codon) -> u8 {
        // TODO drop - at the end (and some *?)
        let &(start, codon) = self.table.get(codon).unwrap_or(&(false, b'-'));
        if start && methionine { b'M' } else { codon }
    }

    /// Print this translation table in a human readable format.
    pub fn print(&self) {
        println!("{}={}", self.name, self.index);
        let mut output: Vec<[u8; 5]> = Vec::new();
        for codon in CodonIterator::new() {
            let (mm, aa) = *self.table.get(&codon).unwrap();
            output.push([aa, if mm { b'M' } else { b'-' }, codon.0, codon.1, codon.2]);
        }
        for (i, name) in ["AAs", "Starts", "Base1", "Base2", "Base3"].into_iter().enumerate() {
            println!("{:<6} = {}", name, output.iter().map(|t5| char::from(t5[i])).collect::<String>());
        }
    }
}

impl str::FromStr for &'static TranslationTable {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        let id = try!(s.parse::<usize>());
        TABLES[id - 1].as_ref().ok_or(ErrorKind::UnknownTable(id).into())
    }
}

/// Represents a DNA codon.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Codon(u8, u8, u8);

impl<'a> From<&'a [u8; 3]> for Codon {
    fn from(bases: &[u8; 3]) -> Self {
        Codon(bases[0], bases[1], bases[2])
    }
}

impl<'a> From<&'a [u8]> for Codon {
    fn from(value: &[u8]) -> Self {
        Codon(value[0], value[1], value[2])
    }
}

struct CodonIterator(usize);

impl CodonIterator {
    fn new() -> Self {
        CodonIterator(0)
    }
}

impl Iterator for CodonIterator {
    type Item = Codon;
    fn next(&mut self) -> Option<Self::Item> {
        if self.0 >= 64 { return None; }
        let next = Codon(BASES[self.0 / 16], BASES[(self.0 / 4) % 4], BASES[self.0 % 4]);
        self.0 += 1;
        Some(next)
    }
}

error_chain! {
    foreign_links {
        ParseTableNumber(std::num::ParseIntError) #[doc = "Indicates failure to parse the table number as a number"];
    }
    errors {
        /// Unknown table number
        UnknownTable(id: usize) {
            description("Table is unknown")
            display("Unknown table: {}", id)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate() {
        assert_eq!(b'L', TABLES[0].as_ref().unwrap().translate(false, &Codon(b'T', b'T', b'G')));
        assert_eq!(b'M', TABLES[0].as_ref().unwrap().translate(true, &Codon(b'T', b'T', b'G')));
    }

    #[test]
    fn test_translation_table_parsing() {
        for (i, mtable) in TABLES.iter().enumerate() {
            let index = format!("{}", i + 1);
            let result = index.parse::<&TranslationTable>();
            match *mtable {
                Some(ref table) => {
                    assert!(result.is_ok());
                    assert_eq!(table, result.unwrap());
                },
                None => {
                    assert!(result.is_err());
                    assert_matches!(*result.unwrap_err().kind(), ErrorKind::UnknownTable(_));
                },
            }
        }
    }

    #[test]
    fn test_codon_from_array() {
        assert_eq!(Codon(b'A', b'A', b'A'), Codon::from(&[b'A', b'A', b'A']));
    }

    #[test]
    fn test_codon_from_slice() {
        assert_eq!(Codon(b'A', b'A', b'A'), Codon::from("AAA".as_bytes()));
    }

    #[test]
    fn test_number_of_codons() {
        assert_eq!(64, CodonIterator::new().count());
    }

}
