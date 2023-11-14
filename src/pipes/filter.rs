//! The [filter](crate::cli::filter::Filter) pipe.

use std::collections::HashSet;

use crate::errors;
use crate::interfaces::fragments::Fragments;

/// An iterator adapter which yields each input [Fragments] with filtered fragments.
pub fn pipe<I: Iterator<Item = errors::Result<Fragments>>>(
    min_length: usize,
    max_length: usize,
    contains: String,
    lacks: String,
    input: I,
) -> Filter<I> {
    Filter {
        min_length,
        max_length,
        contains: contains.chars().collect::<HashSet<char>>(),
        lacks: lacks.chars().collect::<HashSet<char>>(),
        input,
    }
}

/// An iterator over filtered [Fragments].
pub struct Filter<I: Iterator<Item = errors::Result<Fragments>>> {
    min_length: usize,
    max_length: usize,
    contains: HashSet<char>,
    lacks: HashSet<char>,
    input: I,
}

impl<I: Iterator<Item = errors::Result<Fragments>>> Iterator for Filter<I> {
    type Item = errors::Result<Fragments>;

    fn next(&mut self) -> Option<Self::Item> {
        let fragments = self.input.next()?;
        match fragments {
            Err(err) => Some(Err(err)),
            Ok(Fragments { header, peptide, fragments }) => {
                let filtered = fragments.into_iter().filter(|&(from, length)| {
                    if length < self.min_length { return false; }
                    if length > self.max_length { return false; }
                    let set = peptide[from..from+length].chars().collect::<HashSet<char>>();
                    self.contains.intersection(&set).count() == self.contains.len()
                        && self.lacks.intersection(&set).count() == 0
                }).collect();
                Some(Ok(Fragments {
                    header,
                    peptide,
                    fragments: filtered,
                }))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::errors::Error;
    use crate::errors::ErrorKind::Test;

    #[test]
    fn test_empty_input() {
        let input = vec![].into_iter();
        let mut filter = pipe(0, 5, "".to_string(), "".to_string(), input);
        assert_matches!(filter.next(), None);
    }

    #[test]
    fn test_bad_input() {
        let input = vec![Err(Test("wrong input").into())].into_iter();
        let mut filter = pipe(0, 5, "".to_string(), "".to_string(), input);
        assert_matches!(filter.next(), Some(Err(Error(Test("wrong input"), _))));
        assert_matches!(filter.next(), None);
    }

    #[test]
    fn test_filter() {
        let input = vec![
            Ok(Fragments::new("a", "AAAABAAAAB", vec![(0, 1), (0, 4), (0, 5), (0, 10)])),
            Ok(Fragments::new("b", "AAAACAAAAC", vec![(0, 4), (0, 5), (5, 4), (5, 5)])),
            Ok(Fragments::new("c", "BBBBABBBBA", vec![(0, 4), (0, 5), (5, 4), (5, 5)])),
        ].into_iter();
        let mut filter = pipe(2, 5, "A".to_string(), "C".to_string(), input);
        assert_eq!(
            filter.next().unwrap().unwrap(),
            Fragments::new("a", "AAAABAAAAB", vec![(0, 4), (0, 5)])
        );
        assert_eq!(
            filter.next().unwrap().unwrap(),
            Fragments::new("b", "AAAACAAAAC", vec![(0, 4), (5, 4)])
        );
        assert_eq!(
            filter.next().unwrap().unwrap(),
            Fragments::new("c", "BBBBABBBBA", vec![(0, 5), (5, 5)])
        );
        assert_matches!(filter.next(), None);
    }
}
