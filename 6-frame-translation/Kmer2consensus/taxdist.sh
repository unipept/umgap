#/bin/bash

trueTaxID=$1
taxID=$2

trueLineage="$(efetch -db taxonomy -id $trueTaxID -format xml | xtract -pattern Taxon -element Lineage)"
lineage="$(efetch -db taxonomy -id $taxID -format xml | xtract -pattern Taxon -element Lineage)"

java TaxonomicDistance "$trueLineage" "$lineage"
	
