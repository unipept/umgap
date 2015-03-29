#!/bin/bash
awk 'BEGIN {i=0} /^>/ {print i; i = 0; print $0} !/^>/ {i+=1} END {print i}' $1 | tail -n+2
