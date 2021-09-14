//! Some utils.

/// Interleaving iterator.
pub struct Zip<E, I: Iterator<Item = E>> {
    parts: Vec<I>,
}

impl<E, I: Iterator<Item = E>> Zip<E, I> {
    /// Constructor for Zip.
    pub fn new(parts: Vec<I>) -> Self {
        Zip { parts }
    }
}

impl<E, I: Iterator<Item = E>> Iterator for Zip<E, I> {
    type Item = Vec<E>;

    fn next(&mut self) -> Option<Self::Item> {
        self.parts.iter_mut().map(|part| part.next()).collect()
    }
}
