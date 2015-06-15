"""Given a refseq assembly accession ID, return a list of accession ids
from all the genomes in the assembly"""

import sys

from Bio import Entrez

from utils import get_nucleotide_ids


if __name__ == '__main__':
    # Variables
    Entrez.email = "tom.naessens@ugent.be"

    if len(sys.argv) != 2:
        print("Usage: {} [refseq_id]".format(sys.argv[0]))
        exit(1)

    # Get all the proteins
    for nucleotide_id in get_nucleotide_ids(sys.argv[1]):
        print(nucleotide_id)
