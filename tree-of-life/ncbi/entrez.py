import subprocess

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez


# Variables
Entrez.email = "tom.naessens@ugent.be"
PROJECT = 58731

def get_nucleotide_ids(project_id):
    """Get all the nucleotides for a projectid"""
    handle = Entrez.esearch(db="nucleotide", retmax=10, term="{}[BioProject]".format(project_id))
    record = Entrez.read(handle)
    handle.close()

    return record


def get_nucleotide_features(nucleotide_id):
    """Yields all the features of a speficic nucleotide"""
    # Get all we need
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, retmode="text", rettype="gbwithparts")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record.features


def get_nucleotide_translation(nucleotide_id):
    """Gets the translation CDS of a specific nucleotide"""
    for feature in get_nucleotide_features(nucleotide_id):
        if feature.type == 'CDS':
            if 'translation' in feature.qualifiers:
                yield process_translation(feature.qualifiers['translation'][0])


def process_translation(translation):
    result = subprocess.Popen(
        "echo {} | prot2pept | peptfilter | unipept pept2lca".format(translation),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    for pept in result.stdout.readlines():
        print(pept)
    print()


# Parse all the IDs of the nucleotides
for nucleotide_id in get_nucleotide_ids(PROJECT)['IdList']:
    for translation in get_nucleotide_translation(nucleotide_id):
        process_translation(translation)
