import sys

from lca_calculators.tree_based_lca import Tree_LCA_Calculator

class FastaRecord:
    def __init__(self, fasta_id):
        self.fasta_id = fasta_id
        self.peptides = dict()

    def add_peptide(self, peptide, taxon_id):
        if peptide not in self.peptides:
            self.peptides[peptide] = []

        self.peptides[peptide].append(int(taxon_id))


fasta2record = dict()
next(sys.stdin)
for line in sys.stdin:
    fasta_header, peptide, _, taxon_id = line.strip().split(',')

    fasta_id = fasta_header[1:]

    if fasta_id not in fasta2record:
        fasta2record[fasta_id] = FastaRecord(fasta_id)

    fasta2record[fasta_id].add_peptide(peptide, taxon_id)

calculator = Tree_LCA_Calculator()
print("fasta_header,taxon_id,taxon_name,taxon_rank")
for fasta_record in fasta2record.values():

    lcas = [calculator.calc_lca(fasta_record.peptides[peptide]) for peptide in fasta_record.peptides]
    lca = calculator.calc_lca(lcas)

    taxon = calculator.tree.taxons[lca]
    print(">{},{},{},{}".format(fasta_record.fasta_id, taxon.taxon_id, taxon.name, taxon.rank))
