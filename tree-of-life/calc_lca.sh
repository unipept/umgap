#!/bin/bash

cat sample7.txt | unipept pept2prot -s taxon_id | python dosomething.py
