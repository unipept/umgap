import os
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

    unidentified = 0

    FASTA_record = 0
    FASTA_identifier, prev_FASTA_identifier = None, None

    csvcontent = ''
    for line in fileobject:
        line = line.strip()
        if line.startswith('>'):

            FASTA_record += 1
            FASTA_identifier = '_'.join(line.split('_')[:-3])[1:]
            
            if not aggregateDNA or FASTA_identifier != prev_FASTA_identifier:
                
                if csvcontent:
                    try:
                        tree = generateTree(StringIO(csvheader + csvcontent))
                        yield tree
                    except:
                        pass
                else:
                    unidentified += 1
                    
                csvcontent = ''
                
        else:
            
            if line.endswith(',1,root,1,no rank'):
                line += ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,'
                
            csvcontent += '\n' + line
                        
    # return the last tree
    if csvcontent:
        yield generateTree(StringIO(csvheader + csvcontent))
    else:
        unidentified += 1
        
    # return tree that bundles all unassigned DNA/protein reads
    if includeUnidentified:
        tree = TaxonTree()
        lineage = Lineage([Taxon(-1, 'unidentified', 'no rank')])
        tree.addLineage(lineage, unidentified)
        yield tree
    
# merge taxonomic assignments into a single taxonomic tree
handle = open(os.path.join('data', 'sample.prot2pept2lca'), 'r') 
tree = merge_trees(
    generateTrees(
        handle, 
        includeUnidentified=True, 
        aggregateDNA=True
    )
)
handle.close()
#tree.pruneNodes(0.0001)


#tree = merge_trees(generateTrees('test.csv'))
#tree.normalizeChildren()

# focus on a specific node in the tree and rescale the subtree at that node
#tree = tree['Bacteria']
#tree = tree['Eukaryota']
#tree = tree['Archaea']
#tree = tree['Viruses']
#tree.normalizeNode()

# modify HTML visualization template
handle = open(os.path.join('templates', 'spacetree.html'), 'r')
code = handle.read().replace(
    '<JSON_TREE>', 
    tree.export_json(indent=4, showWeight=True)
) 
handle.close()

# output HTML visualization
protein_viz_file = 'sample.treeview.html'
#protein_viz_file = 'test.html'
handle = open(os.path.join('data', protein_viz_file), 'w')
handle.write(code)
handle.close()
