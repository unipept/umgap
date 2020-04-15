
# UMGAP - Unipept Metagenomics Analysis Pipeline

The Unipept Metagenomics Analysis Pipeline can analyse metagenomic samples and
return a frequency table of the taxons it detected for each read. It is based on
the [Unipept Metaproteomics Analysis Pipeline][Unipept]. Both
tools were developed at the Department of Applied Maths, Computer science and
Statistics at [Ghent University].


## Installation & Setup

1. Install [Rust], according to [their installation instructions][rust-install],
   or use your favourite package manager (e.g. `apt install rustc`). The
   pipeline is developed for the latest stable release, but should work on 1.35
   and higher.

2. Clone this repository and go to the repository root.

   ```sh
   git clone https://github.com/unipept/umgap.git
   cd umgap
   ```

3. Compile and install the UMGAP.

   ```sh
   cargo build --release
   cargo install --path .
   ```

   For a multiuser installation, instead of `cargo install`, use `install` to
   place the umgap program and the wrapper script were all users can reach it:

   ```sh
   sudo install target/release/umgap scripts/umgap-analyse.sh /usr/bin
   ```

   `cargo install` will install the `umgap` command to `~/.cargo/bin` by
   default. Please ensure this directory is in your `$PATH`. You can check if
   the installation was succesful by asking for the version:

   ```sh
   umgap -V
   ```

4. (optional) Install [FragGeneScanPlusPlus] to use as gene predictor in the
   pipeline.

5. Run [`scripts/umgap-setup.sh`](scripts/umgap-setup.sh) to interactively
   configure the UMGAP and download the data files required for some steps of
   the pipeline.

   Depending on which type of analysis you are planning, you will need the
   tryptic index file (less powerfull, but runs on any decent laptop) and the
   9-mer index file (uses about 100GB disk space for storage and as much RAM
   during operation. The exact size depends on the version.)

   Run `sudo scripts/umgap-setup.sh -c /etc/umgap -d <datamap>` instead to share
   the datafiles between users. Make sure the `<datamap>` is accessible for the
   end users.

6. (optional) Analyze some test data! Running

   ```sh
   ./scripts/umgap-analyse.sh -1 testdata/A1.fq -2 testdata/A2.fq -t tryptic-sensitivity -o - | tee output.fa
   ```

   should show you a FASTA-like file with a taxon id per header. If you didn't
   download the tryptic index file but the 9-mer index file, use instead:

   ```sh
   ./scripts/umgap-analyse.sh -1 testdata/A1.fq -2 testdata/A2.fq -o - | tee output.fa
   ```

7. (optional) NOT YET INTEGRATED - Visualize some test data! Running

   ```sh
   ./scripts/umgap-visualize.sh output.fa output.html
   ```

   will give you an HTML-file, which will show you a visualization of the test
   data in your favorite browser.

### Updating

A source install can be updated by pulling the repository to get the latest
changes and running in the repository root:

```sh
cargo install --force --path .
```


## Usage

The UMGAP offers individual tools which integrate into a pipeline. Running
`umgap help` will get you to the documentation of each tool, and the short
[metagenomics casestudy] at the Unipept website displays their usage.

This repository also offers 6 preconfigured pipelines
([`scripts/umgap-analyse.sh`](scripts/umgap-analyse.sh))
which should cover most usecases. Running the script without any arguments
should get you started. The preconfigured pipelines are:

- `high-precision`: the default, focusses on high precision with very decent
  sensitivity.
- `max-precision`: focusses on very high precision, at the cost of sensitivity.
- `high-sensitivity`: focusses on high sensitivity, with decent precision.
- `max-sensitivity`: focusses on very high sensitivity, at the cost of
  precision.
- `tryptic-precision`: focusses on high precision, using a much smaller index
  file, which makes it usable on a laptop.
- `tryptic-sensitivity`: focusses on high sensitivity, using a much smaller
  index file, which makes it usable on a laptop.

Another script, [`scripts/umgap-visualize.sh`](scripts/umgap-visualize.sh) will
help you to visualize the output of the pipeline. Again, running the script
without any arguments prints the usage instructions.


## Contributing

Please adhere to the [editorconfig] and [RustFMT] styles specified.


## License

The UMGAP is released under the terms of the MIT License. See the LICENSE file
for more info.


[Unipept]: https://unipept.ugent.be/
[Ghent University]: https://www.ugent.be/
[Rust]: https://www.rust-lang.org/
[metagenomics casestudy]: https://unipept.ugent.be/clidocs/casestudies/metagenomics
[rust-install]: https://www.rust-lang.org/tools/install
[FragGeneScanPlusPlus]: https://github.com/unipept/FragGeneScanPlusPlus
[EditorConfig]: https://editorconfig.org/
[RustFMT]: https://github.com/rust-lang/rustfmt
