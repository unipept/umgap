//! The `umgap prot2kmer2lca` command.

use std::fs;
use std::io;
use std::io::Read;
use std::io::Write;
use std::os::unix::net::UnixListener;
use std::path::PathBuf;
use std::sync::Mutex;

use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::errors;
use crate::io::fasta;

/// The `umgap prot2kmer2lca` command takes one or more peptides as input and outputs the lowest
/// common ancestors of their k-mers.
///
/// The input is given on *standard input* in a FASTA format. Per FASTA header should be a single
/// peptide, which may be hardwrapped with newlines.
///
/// All overlapping k-mers in these peptides (*k* configurable via the `-k` option, and 9 by
/// default) are searched for in the FST index passed as argument. The results are printed out.
///
///     $ cat input.fa
///     >header1
///     DAIGDVAKAYKKAG*S
///     $ umgap prot2kmer2lca -k9 uniprot-2020-04-9mer.index < input.fa
///     >header1
///     571525
///     571525
///     6920
///     6920
///     1
///     6920
///
/// Add the `-o` option to print out 0 for k-mers not found in the index.
///
///     $ umgap prot2kmer2lca -o uniprot-2020-04-9mer.index < input.fa
///     >header1
///     571525
///     571525
///     6920
///     6920
///     1
///     6920
///     0
///     0
///
/// This command also allows an alternative mode of operation. When left on disk, it can take some
/// time for the index to be searched. With the `-m` flag, the complete index will be loaded in
/// memory before operation. This, too, takes some time, but for a single large analysis, this is no
/// hindrance. When processing many short files, the index would need to be loaded again and again.
/// Instead of using this command as part of a pipeline, `... | umgap prot2kmer2lca index | ...`, it
/// can run in a separate (and persistent) process, reusing the same loaded index. Run `umgap
/// prot2kmer2lca -m -s umgap-socket index` somewhere, and when the index is loaded, change your
/// original pipeline(s) to communicate with the socket using OpenBSD's netcat: `... | nc -NU
/// /path/to/umgap-socket | ...`.
#[derive(Debug, StructOpt)]
pub struct ProtToKmerToLca {
    /// The length of the k-mers in the FST
    #[structopt(short = "k", long = "length", default_value = "9")]
    pub length: usize,

    /// Map unknown sequences to 0 instead of ignoring them
    #[structopt(short = "o", long = "one-on-one")]
    pub one_on_one: bool,

    /// An FST to query
    #[structopt(parse(from_os_str))]
    pub fst_file: PathBuf,

    /// Instead of reading from stdin and writing to stdout, create a Unix
    /// socket to communicate with using OpenBSD's netcat (`nc -NU <socket>`).
    /// This is especially useful in combination with the `--in-memory` flag:
    /// you only have to load the FST in memory once, after which you can query
    /// it without having the loading time overhead each time.
    #[structopt(parse(from_os_str), short = "s", long = "socket")]
    pub socket: Option<PathBuf>,

    /// Load FST in memory instead of mmap'ing the file contents. This makes
    /// querying significantly faster, but requires some time to load the FST
    /// into memory.
    #[structopt(short = "m", long = "in-memory")]
    pub fst_in_memory: bool,

    /// How much reads to group together in one chunk, bigger chunks decrease
    /// the overhead caused by multithreading. Because the output order is not
    /// necessarily the same as the input order, having a chunk size which is
    /// a multiple of 12 (all 6 translations multiplied by the two paired-end
    /// reads) will keep sequences of the same reads together.
    #[structopt(short = "c", long = "chunksize", default_value = "240")]
    pub chunk_size: usize,
}

/// Implements the prot2kmer2lca command
pub fn prot2kmer2lca(args: ProtToKmerToLca) -> errors::Result<()> {
    let fst = if args.fst_in_memory {
        let bytes = fs::read(&args.fst_file)?;
        fst::Map::from_bytes(bytes)?
    } else {
        unsafe { fst::Map::from_path(&args.fst_file) }?
    };
    let default = if args.one_on_one { Some(0) } else { None };
    if let Some(socket_addr) = &args.socket {
        let listener = UnixListener::bind(socket_addr)?;
        println!("Socket created, listening for connections.");
        listener
            .incoming()
            .map(|stream| {
                println!("Connection accepted. Processing...");
                let stream = stream?;
                stream_prot2kmer2lca(
                    &stream,
                    &stream,
                    &fst,
                    args.length,
                    args.chunk_size,
                    default,
                )
            })
            .for_each(|result| match result {
                Ok(_) => println!("Connection finished succesfully."),
                Err(e) => println!("Connection died with an error: {}", e),
            });
        Ok(())
    } else {
        stream_prot2kmer2lca(
            io::stdin(),
            io::stdout(),
            &fst,
            args.length,
            args.chunk_size,
            default,
        )
    }
}

fn stream_prot2kmer2lca<R, W>(
    input: R,
    output: W,
    fst: &fst::Map,
    k: usize,
    chunk_size: usize,
    default: Option<u64>,
) -> errors::Result<()>
where
    R: Read + Send,
    W: Write + Send,
{
    let output_mutex = Mutex::new(output);
    fasta::Reader::new(input, true)
        .records()
        .chunked(chunk_size)
        .par_bridge()
        .map(|chunk| {
            let chunk = chunk?;
            let mut chunk_output = String::new();
            for read in chunk {
                // Ignore empty reads and reads with a length smaller than k
                if let Some(prot) = read.sequence.get(0).filter(|p| p.len() >= k) {
                    chunk_output.push_str(&format!(">{}\n", read.header));
                    let mut lcas = (0..(prot.len() - k + 1))
                        .map(|i| &prot[i..i + k])
                        .filter_map(|kmer| fst.get(kmer).map(Some).unwrap_or(default))
                        .map(|lca| lca.to_string())
                        .collect::<Vec<_>>()
                        .join("\n");
                    if !lcas.is_empty() {
                        lcas.push('\n');
                    }
                    chunk_output.push_str(&lcas);
                }
            }
            // TODO: make this the result of the map
            // and print using a Writer
            output_mutex
                .lock()
                .unwrap()
                .write(chunk_output.as_bytes())?;
            Ok(())
        })
        .collect()
}
