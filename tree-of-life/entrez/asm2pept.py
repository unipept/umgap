"""Given a refseq assembly accession ID, return a fasta output of sorted
and filtered tryptic peptides grouped by the protein accession id"""

import subprocess
import sys

from Bio import Entrez

from utils import get_nucleotide_ids, get_nucleotide_features


def get_nucleotide_translation(nucleotide_id):
    """Gets the translation CDS of a specific nucleotide"""

    for features in get_nucleotide_features(nucleotide_id):
        for feature in features:
            if feature.type == 'CDS':
                if 'translation' in feature.qualifiers:
                    yield (feature.qualifiers['protein_id'][0],
                           feature.qualifiers['translation'][0])


def process_translation(translation):
    """Processes a translation by running prot2pept|peptfilter|sort|uniq"""

    result = subprocess.Popen(
        "echo {} | prot2pept | peptfilter | sort | uniq".format(translation),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    for pept in result.stdout.readlines():
        yield pept.strip().decode('utf-8')


def get_all_peptides(assembly_refseq_accession_id):
    """Main function"""

    results = set()

    # Get all the proteins
    for nucleotide_id in get_nucleotide_ids(assembly_refseq_accession_id):
        for insdc, translation in get_nucleotide_translation(nucleotide_id):
            results.add((insdc, translation))


    # Print them out ordered by insdc
    for insdc, prot in sorted(results, key=lambda tup: tup[0]):
        print(">|{}".format(insdc))
        for pept in process_translation(translation):
            print(pept)


if __name__ == '__main__':
    # Variables
    Entrez.email = "tom.naessens@ugent.be"

    if len(sys.argv) != 2:
        print("Usage: {} [refseq_id]".format(sys.argv[0]))
        exit(1)

    get_all_peptides(sys.argv[1])
