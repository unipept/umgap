""" Given a file of pept2prot2filter entries, gets the LCA for every peptide.

Expects a list grouped first by the FASTA header and then by peptide for speed reasons
Not necessarily ordered, just grouped is fine
"""

import sys

from lca_calculators.tree_based_lca import Tree_LCA_Calculator


# Skip the header
next(sys.stdin)

calculator = Tree_LCA_Calculator()
print("fasta_header,taxon_id,taxon_name,taxon_rank")

# Initialize variables
current_fasta_id = ''
current_peptide = ''
taxon_ids = set()
lcas = set()


# Iterate over the lines
for line in sys.stdin:

    fasta_header, peptide, _, taxon_id = line.strip().split(',')
    fasta_id = fasta_header[1:]

    # New peptide: calculate the set of current peptides and add the lca to the fasta_set
    # When there is a new fasta header, we also want to do this
    if (current_peptide and current_peptide != peptide) or (current_fasta_id and current_fasta_id != fasta_id):
        # Calculate the lca of the taxon_ids for the current peptide
        lca = calculator.calc_lca(taxon_ids)

        # Add it to the lcas for the current fasta id
        lcas.add(lca)

        # Empty the taxon_ids set
        taxon_ids = set()

    # New fasta header: calculate the LCA of the peptides so far and print it out
    if current_fasta_id and current_fasta_id != fasta_id:
        # Calculate the lca of the lcas for the current peptide
        lca = calculator.calc_lca(lcas)

        # Get the taxon we're talking about
        taxon = calculator.tree.taxons[lca]
        print(">{},{},{},{}".format(current_fasta_id, taxon.taxon_id, taxon.name, taxon.rank))

        # Empty the sets for the current fasta header and taxon_ids
        lcas = set()

    # Update the current values
    current_fasta_id = fasta_id
    current_peptide = peptide

    # Add the taxon_id to the set
    taxon_ids.add(int(taxon_id))


# And calculate it once again for the last sets
lca = calculator.calc_lca(taxon_ids)
lcas.add(lca)
lca = calculator.calc_lca(lcas)
taxon = calculator.tree.taxons[lca]
print(">{},{},{},{}".format(current_fasta_id, taxon.taxon_id, taxon.name, taxon.rank))

# Finally clean the calculator
calculator.cleanup()
