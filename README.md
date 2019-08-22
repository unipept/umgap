
# UMGAP - Unipept Metagenomics Analysis Pipeline

The Unipept Metagenomics Analysis Pipeline can analyse metagenomic samples and
return a frequency table of the taxons it detected for each read. It is based on
the [Unipept Metaproteomics Analysis Pipeline][Unipept]. Both
tools were developed at the Department of Applied Maths, Computer science and
Statistics at [Ghent University].


## Installation

### From releases

No proper releases are available yet, please install from source.

### From source

1. Install [Rust], according to [their installation instructions][rust-install].

2. Clone this repository and go to the repository root.

   ```sh
   git clone https://github.com/unipept/umgap.git
   cd umgap
   ```

3. Compile and install the UMGAP.

   ```sh
   cargo install --path .
   ```

   This will install the `umgap` command to `~/.cargo/bin` by default. Please
   ensure this directory is in your `$PATH`.

4. (optional) Install [FragGeneScanPlusPlus] to use as gene predictor in the
   pipeline.

#### Updating

A source install can be updated by pulling the repository to get the latest
changes and running in the repository root:

```sh
cargo install --force --path .
```


## Setup

Some parts of the UMGAP requires (large) data files which are not installed
alongside the command. While it is possible to download or even build these
manually, it is advised to use the [setup script] provided in this repository.
The script will interactively download the latest data files and check your
FragGeneScanPlusPlus installation.


## Usage

The UMGAP offers individual tools which integrate into a pipeline. Running
`umgap help` will get you to the documentation of each tool, and the short
[metagenomics casestudy] at the Unipept website displays their usage.

This repository also offers 6 preconfigured pipelines in the [analyse script],
which should cover most usecases. Running the script without any arguments
should get you started.


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
[setup script]: scripts/umgap-setup.sh
[analyse script]: scripts/umgap-analyse.sh
[EditorConfig]: https://editorconfig.org/
[RustFMT]: https://github.com/rust-lang/rustfmt
