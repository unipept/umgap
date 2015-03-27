import subprocess
import urllib.request, urllib.parse, urllib.error
import os
import sys

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez


# Variables
Entrez.email = "tom.naessens@ugent.be"

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

    yield record.features


if len(sys.argv) != 2:
    print("Usage: {} [refseq_id]".format(sys.argv[0]))
    exit(1)

# Get all the proteins
for nucleotide_id in get_nucleotide_ids(sys.argv[1]):
    for features in get_nucleotide_features(nucleotide_id):
        for feature in features:
            if feature.type == 'CDS':
                print(feature.qualifiers['protein_id'][0])
