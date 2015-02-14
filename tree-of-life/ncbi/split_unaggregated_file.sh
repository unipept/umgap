#!/bin/bash
awk '/>\|/ { close("results/fastafiles/"header".fst"); header = $0; sub(/>\|/, "", header); file = "results/fastafiles/"header".fst"; print $0 > file } !/>\|/ { print $0 >> file}'
