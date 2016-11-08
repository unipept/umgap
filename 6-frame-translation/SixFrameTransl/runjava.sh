#!/bin/bash
mvn exec:java -Dexec.mainClass="com.mycompany.sixframetransl.translate" -Dexec.args="sequence.fasta 11 False" 2> /dev/null > output.test
egrep -v '(INFO|WARNING)' output.test > protein.fasta.test
cat protein.fasta.test | prot2pept |peptfilter > peptides.test
