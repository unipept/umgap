
# UMGAP - Unipept Metagenomics Analysis Pipeline

The Unipept Metagenomics Analysis Pipeline can analyse metagenomic samples and
return a frequency table of the taxons it detected for each read. It is based on
the [Unipept Metaproteomics Analysis Pipeline][Unipept]. Both
tools were developed at the Department of Applied Maths, Computer science and
Statistics at [Ghent University](UGent).


## Building

You will need [Rust] installed. The pipeline can then be installed by issuing
the following commands:

```sh
git clone 'git@github.ugent.be:unipept/unipept-metagenomics-scripts.git' umgap
cd umgap
cargo install --path .
```

This will install the `umgap` command to `~/.cargo/bin` by default. Please
ensure this directory is in your `$PATH`.


## Updating the pipeline

Your installed version of the pipeline can be updated by issuing

```sh
cargo install --force --path .
```

## Running the pipeline

See `umgap help` or have a look at the short [metagenomics casestudy] at the
Unipept website.


## License

The UMGAP is released under the terms of the MIT License. See the LICENSE file
for more info.


[Unipept]: https://unipept.ugent.be/
[UGent]: https://www.ugent.be/
[Rust]: https://www.rust-lang.org/
[metagenomics casestudy]: https://unipept.ugent.be/clidocs/casestudies/metagenomics

