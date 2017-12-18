# Test dataset
This subdirectory contains a small dataset to allow tests on the UMGAP pipeline.
The original data comes from [An evaluation of the accuracy and speed of
metagenome analysis tools][metabenchmark] by Lindgreen, Adair & Gardner (2016).
This article is also commonly referred to as the Metabenchmark. The data that is
contained in this directory was randomly sampled (with respect to the different
categories of organisms and their frequencies).

The 2 data files represent paired-end reads. As such, the n'th record in `A1.fq`
will originate from the same organism as the n'th record in `A2.fq`.

## Assigning the taxon ID
Since the original dataset does not contain the taxon ID's of the originally
sampled organisms, we prepended these at the beginning of each record header.
This was done with the help of the `acc2taxonid.py` script. This script takes 2
arguments: the input FASTQ file and a file which maps accession ID's to taxon
ID's. This file can be [downloaded][accession2taxid-file] from the appropriate
NCBI [FTP directory][accestion2taxid-dir].

If a record contains a randomly shuffled nucleotide sequence, the script will
assign it a taxon ID `0`. If the script does not recognize the header in any
way, it will assign a taxon ID `-1`.


[metabenchmark]: http://www.ucbioinformatics.org/metabenchmark.html
[accession2taxid-file]: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
[accestion2taxid-dir]: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
