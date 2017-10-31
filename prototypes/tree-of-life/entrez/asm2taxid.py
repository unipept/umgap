"""Given a refseq assembly accession ID, return the taxon id"""

import sys

from Bio import Entrez

from utils import get_taxon_id


if __name__ == '__main__':
    # Variables
    Entrez.email = "tom.naessens@ugent.be"

    if len(sys.argv) != 2:
        print("Usage: {} [refseq_id]".format(sys.argv[0]))
        exit(1)

    # Get  the taxon id
    print(get_taxon_id(sys.argv[1]))
