#!/bin/bash
for i in 1 2 3 4 5 6 9 10 11
do
	mvn exec:java -Dexec.mainClass="com.mycompany.sixframetransl.translate" -Dexec.args="$1 $i False" 2> /dev/null > output
	egrep -v '(INFO|WARNING)' output > $2.$i.sixframe
	cat $2.$i.sixframe | prot2pept | peptfilter |unipept pept2lca > found_lcas_$2.$i.txt
	mv $2.$i.sixframe results/$3
	mv found_lcas_$2.$i.txt results/$3
done
