"""Utils to query Entrez for stuff"""

import urllib.request

from Bio import SeqIO
from Bio import Entrez


def get_taxon_id(assembly_id):
    """Get the taxon ID for an assembly_id"""
    url = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/" +
           "{}.assembly.txt".format(assembly_id))

    resp = urllib.request.urlopen(url)
    return [line.decode('utf-8').split()[-1]
            for line in resp
            if line.startswith(b'# Taxid')][0]


def get_nucleotide_ids(assembly_id):
    """Get all the nucleotides for a projectid"""

    url = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/" +
           "{}.assembly.txt".format(assembly_id))

    resp = urllib.request.urlopen(url)

    return [line.decode('utf-8').split()[4]
            for line in resp
            if not line.startswith(b'#')]


def get_nucleotide_features(nucleotide_id):
    """Yields all the features of a speficic nucleotide"""

    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id,
                           retmode="text", rettype="gbwithparts")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    yield record.features
