#!/bin/bash

cat data/sample7.txt | unipept pept2prot -s taxon_id | python calc_lca.py
