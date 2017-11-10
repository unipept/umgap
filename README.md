UMGAP - Unipept Metagenomics Analysis pipeline
==============================================
The Unipept Metagenomics Analysis Pipeline can analyse metagenomic samples and
return a frequency table of the taxons it detected for each read. It is based on
the [Unipept Metaproteomics Analysis Pipeline](Unipept). Both
tools were developed at the Department of Applied Maths, Computer science and
Statistics at [Ghent University](UGent).

## Building
To download and build the pipeline, you will need the following programs:

* [Bash] (version 4.0.0 or higher)
* [Meson]
* [GNU coreutils]
* [Rust]
* [Ruby] (and the Unipept gem)

If you have these programs installed, you can build the pipeline by issuing the
following commands:

```
git clone 'git@github.ugent.be:unipept/unipept-metagenomics-scripts.git' umgap
cd umgap
meson build
ninja -C build
```

The pipeline can then be executed with the following command:

```
build/pipeline.sh [OPTIONS] <INPUT FILES>
```

To get a list of options, please run the pipeline with the `--help` option.

## License
The UMGAP is released under the terms of the MIT License. See the LICENSE file
for more info.


[Unipept]: https://unipept.ugent.be/
[UGent]: https://www.ugent.be/
[Bash]: https://www.gnu.org/software/bash/
[Meson]: https://mesonbuild.com/
[GNU coreutils]: https://www.gnu.org/software/coreutils/
[Rust]: https://www.rust-lang.org/
[Ruby]: https://www.ruby-lang.org/
