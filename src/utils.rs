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

/// The bytes backing an FST index.
pub enum FstData {
    /// The whole index, read into memory.
    Memory(Vec<u8>),
    /// The index file, memory mapped.
    Mapped(memmap2::Mmap),
}

impl AsRef<[u8]> for FstData {
    fn as_ref(&self) -> &[u8] {
        match self {
            FstData::Memory(bytes) => bytes,
            FstData::Mapped(mmap) => mmap,
        }
    }
}

/// Opens an FST index, reading it into memory if `in_memory` is set and memory mapping
/// it otherwise.
///
/// fst dropped its own mmap support in 0.4, so this carries the same caveat its
/// `Map::from_path` did: the file must not be modified while the index is in use.
pub fn load_fst<P: AsRef<std::path::Path>>(
    path: P,
    in_memory: bool,
) -> crate::errors::Result<fst::Map<FstData>> {
    let data = if in_memory {
        FstData::Memory(std::fs::read(path)?)
    } else {
        FstData::Mapped(unsafe { memmap2::Mmap::map(&std::fs::File::open(path)?)? })
    };
    Ok(fst::Map::new(data)?)
}
