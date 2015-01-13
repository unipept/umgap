from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez

# Variables
Entrez.email = "tom.naessens@ugent.be"

# Get the project nucleotides
handle = Entrez.esearch(db="nucleotide", retmax=10, term="58731[BioProject]")
record = Entrez.read(handle)
handle.close()

# Parse all the IDs of the nucleotides
for nucleotide_id in record['IdList']:
    # Get all we need
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, retmode="text", rettype="gbwithparts")
    rec = SeqIO.read(handle, "gb")
    handle.close()

    # For every feature of the nucleotide, get the translated sequence
    for feature in rec.features:
        if 'translation' in feature.qualifiers:
            print(feature.qualifiers['translation'][0])
