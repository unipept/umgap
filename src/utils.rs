//! Some utils.
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::os::unix::net::UnixListener;
use std::path::PathBuf;

use errors::Result;
use io::fasta::{transform_records, Record};

/// Interleaving iterator.
pub struct Zip<E, I: Iterator<Item = E>> {
	parts: Vec<I>,
}

impl<E, I: Iterator<Item = E>> Zip<E, I> {
	/// Constructor for Zip.
	pub fn new(parts: Vec<I>) -> Self {
		Zip { parts: parts }
	}
}

impl<E, I: Iterator<Item = E>> Iterator for Zip<E, I> {
	type Item = Vec<E>;

	fn next(&mut self) -> Option<Self::Item> {
		self.parts.iter_mut().map(|part| part.next()).collect()
	}
}

/// Start a daemon by creating a socket at the given socket address, which
/// transforms data parallellised from incoming connections using the
/// `io::fasta::transform_records` function.
pub fn daemonize<F>(socket_addr: &PathBuf, transform: F, chunk_size: usize) -> Result<()>
	where F: Fn(Record) -> Record + Sync {
	let listener = UnixListener::bind(socket_addr)?;
	println!("Socket created, listening for connections.");
	listener.incoming()
	        .map(|stream| {
		        println!("Connection accepted. Processing...");
		        let stream = stream?;
		        transform_records(&stream, &stream, &transform, chunk_size)
		       })
	        .for_each(|result| match result {
		        Ok(_) => println!("Connection finished succesfully."),
	            Err(e) => println!("Connection died with an error: {}", e),
	        });
	unreachable!();
}

/// Set the amount of threads rayon may use.
pub fn set_num_threads(num: usize) -> Result<()> {
	rayon::ThreadPoolBuilder::new().num_threads(num)
	                               .build_global()?;
	Ok(())
}

/// Reads an index file and returns the resulting hashmap. The file is read
/// line by line and each line is split on the given delimiter. The keys will
/// be the first part of the splitted line, the values the second part.
pub fn file2index(file: PathBuf,
                  header: bool,
                  delimiter: &String)
                  -> Result<HashMap<String, String>>
{
	let mut lines = BufReader::new(File::open(file)?).lines();

	// Skip header
	if header {
		lines.next();
	}

	lines.map(|line| {
		     let line = line?;
		     let mut split = line.splitn(2, delimiter);
		     let kmer: String = String::from(split.next().ok_or("Empty line")?);
		     let ids = String::from(split.next().ok_or("No second value")?);
		     Ok((kmer, ids))
		    })
	     .collect()
}

/// Count how much each item in the given iterator occurs, returns a vector
/// with tuples (item, count) sorted descending by count (e.g. most occuring
/// item first).
pub fn item_counts<I, V>(iter: I) -> Vec<(V, usize)>
	where I: Iterator<Item = V>,
	      V: std::cmp::Eq + std::hash::Hash + std::cmp::Ord
{
	// Use a BTreeMap for determinism
	let mut counts: BTreeMap<V, usize> = BTreeMap::new();
	for item in iter {
		let counter = counts.entry(item).or_insert(0);
		*counter += 1;
	}
	let mut counted = counts.into_iter().collect::<Vec<(V, usize)>>();
	counted.sort_unstable_by_key(|(_id, count)| *count);
	counted.reverse();
	return counted;
}
