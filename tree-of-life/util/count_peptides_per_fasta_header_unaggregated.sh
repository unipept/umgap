#!/bin/bash
 tail -n+1 $1 | awk 'BEGIN {i=0} /^>/ {printf("%d\n%s ", i, $0); i = 0;} !/^>/ {i+=1} END {print i}' | tail -n+2 | sort -k1 -t' '
