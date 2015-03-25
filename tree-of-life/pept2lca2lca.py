import sys

from lca_calculators.tree_based_lca import Tree_LCA_Calculator

class FastaRecord:
    def __init__(self, fasta_id):
        self.fasta_id = fasta_id
        self.taxons = []

    def add_taxon(self, taxon):
        self.taxons.append(int(taxon))


fasta2record = dict()
next(sys.stdin)
for line in sys.stdin:
    fasta_header, _, taxon_id, _, _ = line.strip().split(',')

    fasta_id = fasta_header[1:]

    if fasta_id not in fasta2record:
        fasta2record[fasta_id] = FastaRecord(fasta_id)

    fasta2record[fasta_id].add_taxon(taxon_id)

calculator = Tree_LCA_Calculator()
print("fasta_header,taxon_id,taxon_name,taxon_rank")
for fasta_record in fasta2record.values():
    lca = calculator.calc_lca(fasta_record.taxons)

    taxon = calculator.tree.taxons[lca]
    print(">{},{},{},{}".format(fasta_record.fasta_id, taxon.taxon_id, taxon.name, taxon.rank))
