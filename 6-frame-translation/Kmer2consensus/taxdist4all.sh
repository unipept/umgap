#/bin/bash

while read line
do
	header=$(echo $line | cut -d',' -f1 | cut -d'.' -f1)
	taxID=$(echo $line | cut -d',' -f2)
	trueID=$(grep -e "$header.*" Kraken_aanwezige_TaxID.txt |cut -d',' -f3)
	echo "$line,$(./taxdist.sh $taxID $trueID)" >> temp
done < Kraken_HiSeq_consensus.csv
