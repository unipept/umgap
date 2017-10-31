import os, sys
from io import StringIO
from taxonomic_aggregation import TaxonTree, Taxon, Lineage, read_lineages, merge_trees

def generateTree(fileobject, peptideLengthWeighing=True):

    # read all lineages    
    lineages = read_lineages(fileobject)
    
    # compute total length of mapped peptides: longer peptides are less likely
    # to occur by change alone, and should this be given more weight, since it
    # is less likely that they will lead to false identifications; however, on
    # the other hand, the long the peptide, the its likelihood of being properly
    # covered in the sequence databases is lower
    # TODO: this point should be given more attention
    if peptideLengthWeighing:
        peptide_coverage = sum(
            len(peptide) 
            for peptide, _ in lineages
        )

    # merge lineages into a single hierarchy
    tree = TaxonTree()
    for peptide, lineage in lineages: 
        tree.addLineage(
            lineage, 
            len(peptide) / peptide_coverage if peptideLengthWeighing 
            else 1.0 / len(lineages)
        )
        
    # redistribute weights top-down and bottom-up
    tree.redistributeWeightsTopDown()
    tree.percolateWeightsBottomUp()

    # normalize so that root weight equals one (guarantees that all proteins in 
    # the sample have equal weight); this is no longer needed, since lineages
    # are now weighed directly when aggregating them into a single tree
    #tree.normalizeWeights()
    
    return tree

def generateTrees(fileobject, includeUnidentified=True, aggregateDNA=True):
    
    """
    includeUnidentified: indicates whether or not a child should be added to the
        root level, that accumulates all DNA or protein reads for which not a 
        single identification could be made
        
    aggregateDNA: by default, taxa are aggregated for each FASTA record; this
        option allows to do bundle the taxa from successive FASTA records, if 
        they share a common identifier; since there is not agreed upon format
        for the FASTA header line, we have extracted an identifier here that
        merges FASTA records that come from the same DNA read
    """
    
    # defined CSV header (header was removed from the individual FASTA entries
    # to save disk space)
    csvheader = 'sequence,taxon_id,taxon_name,taxon_parent_id,taxon_rank,superkingdom_id,superkingdom_name,kingdom_id,kingdom_name,subkingdom_id,subkingdom_name,superphylum_id,superphylum_name,phylum_id,phylum_name,subphylum_id,subphylum_name,superclass_id,superclass_name,class_id,class_name,subclass_id,subclass_name,infraclass_id,infraclass_name,superorder_id,superorder_name,order_id,order_name,suborder_id,suborder_name,infraorder_id,infraorder_name,parvorder_id,parvorder_name,superfamily_id,superfamily_name,family_id,family_name,subfamily_id,subfamily_name,tribe_id,tribe_name,subtribe_id,subtribe_name,genus_id,genus_name,subgenus_id,subgenus_name,species_group_id,species_group_name,species_subgroup_id,species_subgroup_name,species_id,species_name,subspecies_id,subspecies_name,varietas_id,varietas_name,forma_id,forma_name'

    FASTA_identifier, new_FASTA_identifier = None, None
    aggragated_records = 0
    unidentified_records = 0

    csvcontent = ''
    for line in fileobject:
        line = line.strip()
        if line.startswith('>'):

            new_FASTA_identifier = line.split('_')[0][1:]
            if FASTA_identifier is None:
                FASTA_identifier = new_FASTA_identifier
            
            if not aggregateDNA or new_FASTA_identifier != FASTA_identifier:
                
                aggragated_records += 1

                if csvcontent:
                    yield generateTree(StringIO(csvheader + csvcontent))
                else:
                    unidentified_records += 1
                    
                csvcontent = ''
                FASTA_identifier = new_FASTA_identifier

        else:
            
            if line.endswith(',1,root,1,no rank'):
                line += ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
                
            csvcontent += '\n' + line
                        
    # return the last tree
    aggragated_records += 1
    if csvcontent:
        yield generateTree(StringIO(csvheader + csvcontent))
    else:
        unidentified_records += 1
        
    # return tree that bundles all unassigned DNA/protein reads
    if includeUnidentified:
        tree = TaxonTree()
        lineage = Lineage([Taxon(0, 'unidentified', 'no rank')])
        tree.addLineage(lineage, float(unidentified_records))
        tree.root.setWeight(float(unidentified_records))
        yield tree
        
# merge taxonomic assignments into a single taxonomic tree
if len(sys.argv) >= 2:
    fileobj = open(sys.argv[1], 'r')
else:
    fileobj = sys.stdin

# TODO: in the current implementation, all DNA/protein reads contribute equally
#       to the overall biodiversity estimation; however, some of the reads have
#       a taxonomic aggregation based on a large number of peptide hits, whereas
#       others are based only on a single hit; if possible, there should be an
#       additional weighing for this, given that a proper weighing scheme can be
#       conceived
tree = merge_trees(generateTrees(fileobj))

if len(sys.argv) >= 2:
    fileobj.close()

# export tree in CSV format
print(tree.export_csv(maxdepth=2, minrank='phylum'))
