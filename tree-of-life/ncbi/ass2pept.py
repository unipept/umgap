import subprocess
import urllib.request
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez


# Variables
Entrez.email = "tom.naessens@ugent.be"
ASSEMBLY_ID = "GCF_000015425.1"

def get_nucleotide_ids(assembly_id):
    """Get all the nucleotides for a projectid"""

    url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/{}.assembly.txt".format(assembly_id)

    f = urllib.request.urlretrieve(url, 'assembly.tmp')

    awk_process = subprocess.Popen(
        "cat assembly.tmp | awk '!/^#/ {print $5}'",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    data = awk_process.stdout.readlines()
    data = [dat.decode('utf-8').strip() for dat in data]

    os.remove("assembly.tmp")

    return data


def get_nucleotide_features(nucleotide_id):
    """Yields all the features of a speficic nucleotide"""

    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, retmode="text", rettype="gbwithparts")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record.features


def get_nucleotide_translation(nucleotide_id):
    """Gets the translation CDS of a specific nucleotide"""

    for feature in get_nucleotide_features(nucleotide_id):
        if feature.type == 'CDS':
            if 'translation' in feature.qualifiers:
                yield feature.qualifiers['translation'][0]


def process_translation(i, translation):
    result = subprocess.Popen(
        "echo {} | prot2pept | peptfilter".format(translation),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    print(">|{}".format(i))
    for pept in result.stdout.readlines():
        print(pept.decode('utf-8').strip())

# Get all the proteins
i = 0
for nucleotide_id in get_nucleotide_ids(ASSEMBLY_ID):
    for translation in get_nucleotide_translation(nucleotide_id):
        process_translation(i, translation)
        i+=1
