#!/bin/bash
tail -n+2 $1 | cut -f1 -d',' | sort | uniq -c | sed 's/^ *//' | sort -k2 -t' ' | sed 's/\(.*\) \(.*\)/\2 \1/'
