# Explorative Benchmark

This directory contains scripts and other files related to an
explorative benchmark performed on UMGAP. The `benchmark` script runs
UMGAP with a high number of different configurations on the data sets
used in the paper introducing Kraken.

The `plot-it.py` script generates SVG files which can be imported in
the `all.svg` file (snapping the white rectangle on the topright grid
corner) to create nice graphs.

## Number of tested configurations

```python
frag = 3
agg = 5
minfreq = 5
kmer = 1
maxgap = 5
minseed = 3
print("6ft", frag * agg * minfreq * (kmer + maxgap * minseed))

minlength = 6
maxlength = 6
print("tryp", frag * agg * minfreq * minlength * maxlength)
```
