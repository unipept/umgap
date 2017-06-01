# Thesis

Repository for all scripts used in my master's thesis

### Adler_scripts
	Contains all bash-scripts used on Adler
	* [kmer.sh](Adler_scripts/kmer.sh): splits a six-frame translation into identified 9-mers
	* [kraken_analysis](Adler_scripts/kraken_analysis): used to make the csv-files used in the R-scripts based on the results obtained on the Kraken-dataset
	* [lineages.py](Adler_scripts/lineages.py): used to make the csv-files used in the R-scripts based on the results obtained on the Metabenchmark-dataset (by F. Van der Jeugt)
	* [score_peptides](Adler_scripts/score_peptides): gives the tryptic peptides scores before passing them on to taxa2agg (works with the new pept2lca)
	* [seedextend](Adler_scripts/seedextend): used to run the seed-extend method
	* [trypept](Adler_scripts/trypept_analysis): contains all comands needed to analyse DNA-reads by splitting them into tryptic peptides	

### Kmer2consensus
	Contains all R-scripts used to analyse results of the tests on the Kraken and Metabenchmarkdataset

### ScoreReads
	Contains the Java code for the plots of the sixframe translations, split into 9mers or tryptic peptides

### Scripts
	Contains all bash-scripts used on my own computer
	* [DrawTikzPlots.sh](Scripts/DrawTikzPlots.sh): Bundels the nessecary preparations and the Java-code to make a visualisation out of a file containing the six-frame translation and a file containing the lca's
	* [getProteins](Scripts/getProteins): retreives all protein positions based on a GenBank-id
	* [score_peptides](Scripts/score_peptides): older version of the script that gives scores to tryptic peptides (works with output from unipept pept2lca)
	* [taxonInfo](Scripts/taxonInfo): retrieves taxon name and rank based on a taxon id
	* [translate_sixframe.sh](Scripts/translate_sixframe.sh): determines the six-frame translation and the lca's of the peptides in that translation

### SeedExtend
	Contains the Java code for the seed-extend method

### Simulaties
	Contains the scripts for the simulations made to determine the score for tryptic peptides

### SixFrameTransl
	Contains the Java code for the six-frame translation

### drawLca4FGS
	Contains the Java code to generate extra code for the plots of the six-frame translation so it can be compared with FragGeneScan output
	
