//! Translates a DNA read into a peptide, using the given frames.

use dna::{Read, Strand, Frame};
use dna::translation::TranslationTable;
use tools::{Record, Peptide};

/// Translates a DNA read into a peptide, using the given frames.
pub fn translate<'a>(table: &'a TranslationTable,
                     frames: &'a Vec<Frame>,
                     methionine: bool,
                     append_name: bool,
                     records: impl IntoIterator<Item = Record<Read>> + 'a)
                     -> impl Iterator<Item = Record<Peptide>> + 'a
{
	records.into_iter()
		.map(|r| (r.header, Strand::from(&r.content)))
		.map(|(h, s)| {
			let r = s.reversed();
			(h, s, r)
		})
		.flat_map(move |(header, forward, reverse)| frames.iter().map(move |frame| Record {
			header: if !append_name { header.clone() } else { header.clone() + "|" + &frame.to_string() },
			content: String::from_utf8(table.translate_frame(methionine, match frame {
				Frame::Forward1 => forward.frame(1),
				Frame::Forward2 => forward.frame(2),
				Frame::Forward3 => forward.frame(3),
				Frame::Reverse1 => reverse.frame(1),
				Frame::Reverse2 => reverse.frame(2),
				Frame::Reverse3 => reverse.frame(3),
			})).unwrap()
		}))
}
