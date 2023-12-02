//! The [bestof](crate::cli::bestof::BestOf) pipe.

use crate::errors;
use crate::interfaces::fragment_ids::FragmentIds;

/// An iterator adapter which yields the best [FragmentIds]
/// of each chunk of `frames` FragmentIds in the input
/// iterator. See also: [bestof].
pub fn pipe<I: Iterator<Item = errors::Result<FragmentIds>>>(
    frames: usize,
    input: I
) -> BestOf<I> {
    BestOf { frames, input }
}

/// An iterator over the best [FragmentIds] found in each chunk of [FragmentIds], made by the [pipe]
/// function.
pub struct BestOf<I: Iterator<Item = errors::Result<FragmentIds>>> {
    frames: usize,
    input: I
}

impl<I: Iterator<Item = errors::Result<FragmentIds>>> Iterator for BestOf<I> {
    type Item = Result<FragmentIds>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk: Vec<FragmentIds> = Vec::with_capacity(self.frames);
        for _c in 0..self.frames {
            match self.input.next() {
                Some(Err(err)) => return Some(Err(err).chain_err(|| "unabled to read one of the frames")),
                Some(Ok(item)) => chunk.push(item),
                None => {},
            }
        }
        if chunk.len() == self.frames {
            Some(Ok(bestof(chunk)))
        } else if chunk.is_empty() {
            None
        } else {
            Some(Err(ErrorKind::MissingFrames(chunk.len(), self.frames).into()))
        }
    }
}

/// Given a chunk of [FragmentIds], return the [FragmentIds] with
/// the most [TaxonId](crate::taxon::TaxonId)s that aren't 0 (no
/// identification) or 1 (root identification).
pub fn bestof(chunk: Vec<FragmentIds>) -> FragmentIds {
    chunk
        .into_iter()
        .max_by_key(|fragments| {
            fragments
                .taxa
                .iter()
                .filter(|&s| *s != 0 && *s != 1)
                .count()
        })
        .unwrap()
}

error_chain! {
    errors {
        /// Missing frames
        MissingFrames(expected: usize, found: usize) {
            description("Missing frames")
            display("Missing frames: expected {}, found {}", expected, found)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bad_input() {
        let mut input = vec![Err(errors::ErrorKind::Test("wrong input").into())].into_iter();
        let mut bestof = pipe(1, &mut input);
        assert_matches!(bestof.next(), Some(Err(Error(ErrorKind::Msg(_), _))));
        assert_matches!(bestof.next(), None);
    }

    #[test]
    fn test_frame_shortage() {
        let mut input = vec![Ok(FragmentIds::new("any|1", vec![0]))].into_iter();
        let mut bestof = pipe(2, &mut input);
        assert_matches!(
            bestof.next(),
            Some(Err(Error(ErrorKind::MissingFrames(1, 2), _)))
        );
        assert_matches!(bestof.next(), None);
    }

    #[test]
    fn test_some_better_than_none() {
        let mut input = vec![
            Ok(FragmentIds::new("any|1", vec![0])),
            Ok(FragmentIds::new("any|2", vec![5])),
        ]
        .into_iter();
        let mut bestof = pipe(2, &mut input);
        let best = bestof.next();
        assert_matches!(best, Some(Ok(_)));
        assert_eq!(best.unwrap().unwrap(), FragmentIds::new("any|2", vec![5]));
        assert_matches!(bestof.next(), None);
    }

    #[test]
    fn test_many_better_than_some() {
        let mut input = vec![
            Ok(FragmentIds::new("any|1", vec![2])),
            Ok(FragmentIds::new("any|2", vec![5, 3, 2, 1])),
        ]
        .into_iter();
        let mut bestof = pipe(2, &mut input);
        let best = bestof.next();
        assert_matches!(best, Some(Ok(_)));
        assert_eq!(
            best.unwrap().unwrap(),
            FragmentIds::new("any|2", vec![5, 3, 2, 1])
        );
        assert_matches!(bestof.next(), None);
    }

    #[test]
    fn test_multiple_chunks() {
        let mut input = vec![
            Ok(FragmentIds::new("any|1", vec![2])),
            Ok(FragmentIds::new("any|2", vec![5, 3, 2, 1])),
            Ok(FragmentIds::new("any|1", vec![2])),
            Ok(FragmentIds::new("any|2", vec![1])),
        ]
        .into_iter();
        let mut bestof = pipe(2, &mut input);
        assert_eq!(
            bestof.next().unwrap().unwrap(),
            FragmentIds::new("any|2", vec![5, 3, 2, 1])
        );
        assert_eq!(
            bestof.next().unwrap().unwrap(),
            FragmentIds::new("any|1", vec![2])
        );
        assert_matches!(bestof.next(), None);
    }
}
