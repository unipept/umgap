//! Compact representation of amino acids

use self::ErrorKind::{InvalidCompactRepresentation, UnknownAminoAcid};
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use num_traits::identities::{One, Zero};

const AA_COUNT: usize = 28;
const COMPACT_BITS: usize = 5;

lazy_static! {
	// Extended IUPAC protein sequence notation
	static ref NUM_AA: [u8; AA_COUNT] = [
		b'-', // 00000 (To prevent AA's with all zeroes)
		b'*', // 00001
		b'A', // 00010
		b'B', // 00011
		b'C', // 00100
		b'D', // 00101
		b'E', // 00110
		b'F', // 00111
		b'G', // 01000
		b'H', // 01001
		b'I', // 01010
		b'J', // 01011
		b'K', // 01100
		b'L', // 01101
		b'M', // 01111
		b'N', // 10000
		b'O', // 10001
		b'P', // 10010
		b'Q', // 10011
		b'R', // 10100
		b'S', // 10101
		b'T', // 10110
		b'U', // 10111
		b'V', // 11000
		b'W', // 11001
		b'X', // 11010
		b'Y', // 11011
		b'Z', // 11100
	];
	static ref AA_NUM: [Option<u8>; 256] = {
		let mut array = [None; 256];
		for i in 0..NUM_AA.len() {
			array[NUM_AA[i] as usize] = Some(i as u8);
		}
		array
	};
}

/// Convert an amino acid sequence to a more compact representation.
pub fn to_compact_1(aa_seq: &[u8]) -> Result<Vec<u8>> {
	let mut compact: BigUint = Zero::zero();
	let mut multiplier: BigUint = One::one();
	let len = aa_seq.len();
	for i in 0..len {
		let aa = aa_seq[len - i - 1];
		let num: u8 = AA_NUM.get(aa as usize)
		                    .unwrap_or(&None)
		                    .ok_or(UnknownAminoAcid(aa))?;
		compact += num * &multiplier;
		multiplier *= AA_COUNT;
	}
	Ok(compact.to_bytes_be())
}

/// Convert an amino acid sequence to a more compact representation.
pub fn to_compact(aa_seq: &[u8]) -> Result<Vec<u8>> {
	let mut current_bit = 0;
	let mut byte: u8 = 0;
	let mut bytes: Vec<u8> = Vec::new();
	for &aa in aa_seq {
		let aa: u8 = AA_NUM.get(aa as usize)
		                   .unwrap_or(&None)
		                   .ok_or(UnknownAminoAcid(aa))?;

		// Lower byte
		byte |= aa << current_bit;

		current_bit += COMPACT_BITS;

		// If we 'overflowed':
		if current_bit >= 8 {
			current_bit -= 8;
			bytes.push(byte);
			byte = aa >> (COMPACT_BITS - current_bit);
		}
	}
	bytes.push(byte);
	Ok(bytes)
}


/// Convert a compact amino acid representation to normal.
pub fn from_compact_1(compact_vec: Vec<u8>) -> Result<Vec<u8>> {
	let mut compact = BigUint::from_bytes_be(&compact_vec);
	let mut aa_seq = Vec::new();
	while !compact.is_zero() {
		let num = &compact % AA_COUNT;
		compact /= AA_COUNT;
		aa_seq.push(NUM_AA[num.to_usize().ok_or("Something went horribly wrong")?]);
	}
	aa_seq.reverse();
	Ok(aa_seq)
}

/// Convert a compact amino acid representation to normal.
pub fn from_compact(compact_vec: Vec<u8>) -> Result<Vec<u8>> {
	let mut aa_count = (compact_vec.len() * 8) / 5;
	if let Some(last_byte) = compact_vec.last() {
		if (last_byte.leading_zeros() as usize) > COMPACT_BITS {
			aa_count -= 1;
		}
	}
	let mut aa_seq = Vec::new();
	let mut current_bit = 0;
	let mut compact = compact_vec.iter();
	if let Some(byte) = compact.next() {
		let mut byte = byte;
		while aa_seq.len() < aa_count {
			// Lower byte
			let mut aa = (byte >> current_bit) & 0b11111;
			current_bit += COMPACT_BITS;

			// If we 'overflowed':
			if current_bit >= 8 {
				current_bit -= 8;
				byte = compact.next()
				              .ok_or(InvalidCompactRepresentation(compact_vec.clone()))?;
				aa |= (byte << (COMPACT_BITS - current_bit)) & 0b11111;
			}
			aa_seq.push(NUM_AA[aa as usize]);
		}
	}
	Ok(aa_seq)
}

error_chain! {
	errors {
		/// Unknown Amino Acid
		UnknownAminoAcid(aa: u8) {
			description("Encountered an unknown amino acid character")
			display("Unknown amino acid: {}", *aa as char)
		}
		/// Invalid compact representation of an amino acid
		InvalidCompactRepresentation(compact: Vec<u8>) {
			description("Encountered an invalid compact representation while trying to convert to normal")
			display("Invalid representation: {:?}", compact)
		}
	}
}

#[cfg(test)]
mod test {
	use super::*;

	const AA_SEQ: &'static str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

	#[test]
	fn test_convert_to() -> () {
		let aa_bytes = Vec::from(AA_SEQ.as_bytes());

		let compact_bytes = to_compact(&aa_bytes);
		assert!(compact_bytes.is_ok());
	}

	#[test]
	fn test_conversion() -> () {
		let aa_bytes = Vec::from(AA_SEQ.as_bytes());

		let compact_bytes = to_compact(&aa_bytes);
		assert!(compact_bytes.is_ok());

		let compact_bytes = compact_bytes.unwrap();
		assert!(compact_bytes.len() < AA_SEQ.len());

		let aa_converted = from_compact(compact_bytes);
		assert!(aa_converted.is_ok());
		let aa_converted = String::from_utf8(aa_converted.unwrap()).unwrap();

		assert_eq!(AA_SEQ, aa_converted);
	}

	#[test]
	fn test_conversion_stars() -> () {
		let aa = "***";
		let aa_bytes = Vec::from(aa.as_bytes());

		let compact_bytes = to_compact(&aa_bytes);
		assert!(compact_bytes.is_ok());

		let compact_bytes = compact_bytes.unwrap();
		assert!(compact_bytes.len() < aa.len());

		let aa_converted = from_compact(compact_bytes);
		assert!(aa_converted.is_ok());
		let aa_converted = String::from_utf8(aa_converted.unwrap()).unwrap();

		assert_eq!(aa, aa_converted);
	}

	#[test]
	fn test_ordening_same_length() -> () {
		let test_seqs: [&'static str; 7] = ["***", "A*A", "AAA", "KAA", "KYA", "YAY", "YYY"];

		let converted = test_seqs.iter()
		                         .map(|seq| to_compact(seq.as_bytes()))
		                         .collect::<Result<Vec<Vec<u8>>>>()
		                         .unwrap();

		for i in 0..(converted.len() - 1) {
			assert!(converted[i] < converted[i + 1]);
		}
	}

}
