#!/usr/bin/env python3
#
# Takes 2 arguments:
# (1) An accession-to-taxon mapping file from the NCBI FTP server
#       (see "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/")
# (2) a FASTQ mapping file from the Metabenchmark
#
# The script returns the FASTQ file, with the header modified in such a way that it now contains
# the taxon id.
#
# NOTE: if you get errors such as "Segmentation Fault" / "Bus Error" / "Memory Error",
#       there's a big chance that you need more memory (as loading the full mapping
#       file takes in quite a lot of space)

import pandas as pd
import re
import sys

# For some god-awful reason the authors of the metabenchmark though it was a good idea to use data
# from the European Neucleotide Archive as well. Since mapping to a proper taxon id isn't always
# possible (or it needs elaborate parsing, and internet access) we'll just hardcode them here.
ENA_mapping = [
    ('AL450380', 272631), # Mycobacterium leprae TN
    ('AL645882', 1902),   # Streptomyces coelicolor
    ('CM000636', 1768),   # Mycobacterium kansasii
    ('CM000789', 1773),   # Mycobacterium tuberculosis
    ('BX548020', 84588),  # Synechococcus sp. WH 8102
    ('CM000748', 527025), # Bacillus thuringiensis serovar thuringiensis
    ('CM000833', 320371), # Burkholderia pseudomallei 1710a
    ('CM000855', 683083), # Campylobacter jejuni subsp. jejuni 414
    ('BX119912', 243090), # Rhodopirellula baltica SH 1
    ('BX470248', 520),    # Bordetella pertussis
    ('BX470250', 518),    # Bordetella bronchiseptica
    ('BX571963', 258594), # Rhodopseudomonas palustris CGA009
    ('BX842601', 50701),  # Bdellovibrio bacteriovorus HD100
    ('CR354532', 74109),  # Photobacterium profundum
    ('AJ235269', 272947), # Rickettsia prowazekii str. Madrid E
    ('AL732656', 211110), # Streptococcus agalactiae NEM316
    ('BX548020', 84588),  # Synechococcus sp. WH 8102
    ('BX571656', 273121), # Wolinella succinogenes DSM 1740
    ('CM000438', 271848), # Burkholderia thailandensis E264
    ('CM000488', 535026), # Bacillus subtilis subsp. subtilis str. NCIB 3610
    ('CM000657', 479833), # Clostridioides difficile QCD-97b34
    ('CM000715', 526980), # Bacillus cereus ATCC 10876
    ('CM000724', 526975), # Bacillus cereus BDRD-ST26
    ('CM000731', 526984), # Bacillus cereus Rock3-29
    ('CM000750', 527027), # Bacillus thuringiensis serovar pakistani str. T13001
    ('CM000754', 527032), # Bacillus thuringiensis serovar andalousiensis BGSC 4AW1
    ('CT009589', 1718),   # Corynebacterium glutamicum
]

def main(args):
    # Read in the mapping table
    acc2taxid = pd.read_csv(args[2], sep='\t', usecols=['accession','taxid'])
    acc2taxid.set_index('accession', inplace=True)

    def get_taxon_id(accession):
        try:
            return acc2taxid.loc[accession]['taxid']
        except (KeyError) as e:
            return -1

    with open(args[1]) as fastq_file:
        for i, line in enumerate(fastq_file):
            # Leave the non-header lines alone
            if i % 4 != 0:
                print(line, end='')
                continue

            # First, take on the edge cases:
            # (Use the second Excel appendix of the Metabenchmark to get an idea of the data)

            # (1) Randomly shuffled records
            if 'random' in line.lower():
                taxon_id = 0

            # (2) The records that were simulated with the Rose simulator.
            #     The paper mentions these can all be traced back to the same accession number
            elif '_divergence_' in line:
                taxon_id = get_taxon_id('AE016823')

            # (3) This is the fun part: some records contain barely any information at all, except
            #     for sometimes a (shortened and unusable for search) organism name and possibly a
            #     chromosome number. This is where you'll need the appendix...
            elif 'Arabidopsis' in line:
                taxon_id = 3702  # Arabidopsis thaliana: a thale cress
            elif 'chr22-_Eukaryotes' in line:
                taxon_id = 9606  # Homo Sapiens
            elif 'chr28-_Eukaryotes' in line:
                taxon_id = 9031  # Gallus gallus: a chicken
            elif 'chr2-_Eukaryotes' in line:
                taxon_id = 28377 # Small caution here: we don't have an organism name, only that it
                                 # should be a lizard. Putting the relevant sequences through BLAST
                                 # gives us the "Anolis carolinensis", i.e. the Carolina anole.
            elif 'Pan_troglo' in line:
                taxon_id = 9598  # Pan troglodytes: a chimpansee
            elif ('JH584390' in line) or ('JH584391' in line):
                taxon_id = 8478  # Chrysemys picta bellii: a western painted turtle
            elif 'Haliaeetus_leucocephalus' in line:
                taxon_id = 52644 # Haliaeetus Leucocephalus: a bald eagle

            # (4) Some records have accessions that belong to ENA instead of NCBI.
            #     Since they don't all map cleany to a single taxon ID, we need manual lookups again.
            elif line.startswith('@ENA|'):
                found = False
                for accession, taxid in ENA_mapping:
                    if line[5:13] == accession:
                        taxon_id = taxid
                        found = True
                if not found:
                    taxon_id = -1   # This shouldn't happen


            # General case: strip the accession number (which comes right after the '@')
            else:
                # One last edge case: some records start with "tpg|" (which can be stripped)
                try:
                    accession_nr = re.search(r'^@(tpg\|)?([A-Z0-9]*)', line).group(2)
                    taxon_id = get_taxon_id(accession_nr)
                except AttributeError:
                    taxon_id = -1   # This shouldn't happen

            print('@{}###{}'.format(taxon_id, line[1:]), end='')

if __name__ == '__main__':
    main(sys.argv)
