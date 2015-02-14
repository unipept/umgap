import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez


# Variables
Entrez.email = "tom.naessens@ugent.be"
PROJECT = 58731

def get_protein_ids(project_id):
    """Get all the proteins for a projectid"""

    handle = Entrez.esearch(db="protein", term="{}[BioProject]".format(project_id), retmax=10000)
    record = Entrez.read(handle)
    handle.close()

    return record

# Get all the proteins
handle = Entrez.efetch(db="protein", id=','.join(get_protein_ids(PROJECT)['IdList']), retmode='xml')
records = Entrez.read(handle, validate=False)
handle.close()

for record in records:
    print(record['GBSeq_accession-version'])
