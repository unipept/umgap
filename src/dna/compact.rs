//! Compact representation of amino acids

use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use num_traits::identities::{One, Zero};

const AA_COUNT: usize = 22;

lazy_static! {
	static ref NUM_AA: [u8; AA_COUNT] = [
	                                     b'-', b'*', b'A', b'C', b'D', b'E', b'F', b'G', b'H',
	                                     b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S',
	                                     b'T', b'V', b'W', b'Y',
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
pub fn to_compact(aa_seq: &[u8]) -> Result<Vec<u8>> {
	let mut compact: BigUint = Zero::zero();
	let mut multiplier: BigUint = One::one();
	let len = aa_seq.len();
	for i in 0..len {
		let aa = aa_seq[len - i - 1];
		let num: u8 = AA_NUM.get(aa as usize)
		                    .unwrap_or(&None)
		                    .ok_or(ErrorKind::UnknownAminoAcid(aa))?;
		compact += num * &multiplier;
		multiplier *= AA_COUNT;
	}
	Ok(compact.to_bytes_be())
}


/// Convert a compact amino acid representation to normal.
pub fn from_compact(compact_vec: Vec<u8>) -> Result<Vec<u8>> {
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

error_chain! {
	errors {
		/// Unknown Amino Acid
		UnknownAminoAcid(aa: u8) {
			description("Encountered an unknown amino acid character")
			display("Unknown amino acid: {}", *aa as char)
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
