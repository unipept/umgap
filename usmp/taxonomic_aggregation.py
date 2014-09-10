def row2lineage(header, row):
    
    ids = [
        int(value) if value else None
        for name, value in zip(header, row)
        if (
            name.endswith('_id') 
            and not name.startswith('taxon_')
        )
    ]
    
    names = [
        value
        for name, value in zip(header, row)
        if (
            name.endswith('_name') 
            and not name.startswith('taxon_')
        )
    ]
    
    ranks = [
        ' '.join(name.replace('__', '_').split('_')[:-1])
        for name, _ in zip(header, row)
        if (
            name.endswith('_name') 
            and not name.startswith('taxon_')
        )
    ]
    
    return (
        row[0],
        Lineage(
            Taxon(*args) if args[0] is not None else None for args in zip(ids, names, ranks)
        )
    )
    
def read_lineages(fileobject):
    
    """
    >>> import os
    >>> lineages = read_lineages(open(os.path.join('data', 'id_protein10.csv'), 'r'))
    >>> lineages
    [('NLLETGNMGR', Taxon(976, 'Bacteroidetes', 'phylum')), ('VVELK', Taxon(1, 'cellular organism', 'root')), ('NGEYIPSFISIDK', Taxon(817, 'Bacteroides fragilis', 'species')), ('LTNEVVAMK', Taxon(171549, 'Bacteroidales', 'order')), ('AENAFIPR', Taxon(817, 'Bacteroides fragilis', 'species')), ('NEGIGK', Taxon(1, 'cellular organism', 'root'))]
    >>> set(lineage for peptide, lineage in lineages)
    {Taxon(976, 'Bacteroidetes', 'phylum'), Taxon(1, 'cellular organism', 'root'), Taxon(171549, 'Bacteroidales', 'order'), Taxon(817, 'Bacteroides fragilis', 'species')}
    """

    import csv    
    csvreader = csv.reader(fileobject, delimiter=',')
    header = next(csvreader)
    return [row2lineage(header, row) for row in csvreader]
    
class Taxon:
    
    """
    >>> taxon = Taxon(976, 'Bacteroidetes', 'phylum')
    >>> print(taxon)
    Bacteroidetes (976, phylum)
    >>> taxon
    Taxon(976, 'Bacteroidetes', 'phylum')
    >>> hash(taxon)
    976
    """
    
    def __init__(self, id, name, rank):
        
        self.id = id
        self.name = name
        self.rank = rank
        
    def getId(self):
        
        return self.id
        
    def getName(self):
        
        return self.name
        
    def getRank(self):
        
        return self.rank
        
    def __str__(self):
        
        return '{} ({}, {})'.format(self.name, self.id, self.rank)

    def __repr__(self):
        
        return 'Taxon({}, {!r}, {!r})'.format(self.id, self.name, self.rank)
    
    def __eq__(self, other):
        
        return self.id == other.id
    
    def __lt__(self, other):
        
        return self.name < other.name
    
    def __le__(self, other):
        
        return self.name <= other.name
    
    def __hash__(self):
        
        return self.id

class Lineage:
    
    def __init__(self, taxa):
        
        self.taxa = list(taxa)
        
    def getTaxa(self):
        
        return self.taxa
        
    def leaf(self):
        
        for taxon in reversed(self.taxa):
            if taxon is not None:
                return taxon
        return Taxon(1, 'cellular organism', 'root')
        
    def __repr__(self):
        
        return repr(self.leaf())
        
    def __eq__(self, other):
        
        return all(
            (
                (taxon1 is None and taxon2 is None)
                or
                (
                    taxon1 is not None and 
                    taxon2 is not None and
                    taxon1 == taxon2
                )
            )
            for taxon1, taxon2 in zip(self.taxa, other.taxa)
        )
        
    def __hash__(self):
        
        return hash(self.leaf())
    
class TaxonTreeNode:
    
    def __init__(self, taxon, weight=0.0):
        
        # set taxon corresponding to current node
        self.taxon = taxon
        
        # set list of child nodes for current node
        self.children = set()
            
        # set weight of current node
        self.weight = weight
        
    def getTaxon(self):
        
        return self.taxon
    
    def getChildren(self):
        
        return self.children
    
    def getWeight(self):
        
        return self.weight
    
    def setWeight(self, weight=0.0):
        
        self.weight = weight
        
    def addChild(self, taxon, weight=0.0):
        
        node = TaxonTreeNode(taxon, weight)
        self.getChildren().add(node)
        return node
    
    def toText(self, branches=None, lastChild=False):
        
        if branches is None:
            branches = []
        
        # show node with appropriate indentation
        indent = ''.join(
            {True:'|   ', False:'    '}[branch] 
            for branch in branches
        )
        text = '{}{}-- {} (weight={})\n'.format(
            indent, 
            {False:'|', True:'`'}[lastChild], 
            self.getTaxon(), 
            self.getWeight()
        )
        
        # adjust branches for child nodes
        if lastChild:
            branches = branches + [False]
        else:
            branches = branches + [True]

        # show child nodes with appropriate indentation
        children = sorted(self.getChildren())
        for node in children[:-1]:
            text += node.toText(branches, lastChild=False)
        if children:
            text += children[-1].toText(branches, lastChild=True)
        return text
            
    def __str__(self):
        
        return self.toText(lastChild=True).rstrip('\n')
    
    def __repr__(self):
        
        return '{} [weight={}]'.format(self.getTaxon(), self.getWeight())
    
    def __eq__(self, other):
        
        return self.getTaxon() == other.getTaxon()
    
    def __hash__(self):
        
        return hash(self.getTaxon())
    
    def __lt__(self, other):
        
        return self.getTaxon() < other.getTaxon()
    
    def nodes(self, includeNode=True):
        
        # return the current node
        if includeNode:
            yield self
            
        # return child nodes and all their subnodes
        for child in self.getChildren():
            for node in child.nodes():
                yield node
        
    def isLeaf(self):
        
        return not bool(self.getChildren())
    
    def leafs(self, includeNode=True):
        
        for node in self.nodes(includeNode):
            if node.isLeaf():
                yield node
                
    def redistributeWeightsTopDown(self):
        
        if self.isLeaf():
            
            # weights of leaf nodes remain as they are and the weight of the 
            # leaf nodes is returned by the function
            return self.getWeight()
        
        else:
            
            # start by redistributing weights of all child nodes (depth-first)
            # and along the way compute the total weight of the subtree, not
            # taking into account the weight of the node itself
            subtree_weight = sum(
                child.redistributeWeightsTopDown() 
                for child in self.getChildren()
            )
            
            # distribute the weight of the node over its leafs; this is not 
            # done evenly, but taking into account the weights of the leafs; 
            # note that at this stage, only the node itself and its leafs have
            # non-zero weights; the internal nodes have already been processed,
            # which has set their weight to zero
            node_weight = self.getWeight()
            for leaf in self.leafs():
                leaf_weight = leaf.getWeight()
                leaf.setWeight(leaf_weight * (1 + node_weight / subtree_weight))

            # weight of node has been distributed over its children
            self.setWeight(0.0)
            
            # return the total weight of the subtree (including the one of the
            # node that has been redistributed); the total weight of the tree
            # remains the same before and after redistribution
            return node_weight + subtree_weight
        
    def getTotalWeight(self, includeNode=True, leafsOnly=False):
        
        return sum(
            node.getWeight() 
            for node in 
            (self.leafs(includeNode) if leafsOnly else self.nodes(includeNode))
        )
        
    def percolateWeightsBottomUp(self):
        
        if self.isLeaf():
            # weights of leaf nodes remain as they are and the weight of the 
            # leaf nodes is returned by the function
            return self.getWeight()
        else:
            # propagation of weights is done first for all child nodes and then
            # the sum of the weights of the child nodes is taken as the new
            # weight of the internal node in the tree
            weight = self.getWeight() + sum(
                child.percolateWeightsBottomUp() for child in self.getChildren()
            )
            self.setWeight(weight)
            return weight
        
    def generate_label(
        node, 
        showWeight=False,
        label_weight_factor=100.0,
        label_weight_formatting='{:.2f}%'
    ):
        
        taxon = node.getTaxon()
        
        label = taxon.getName().split()
        if (
            len(label) > 1 and 
            len(label[0]) > 3 and
            label[0][0].isupper() and
            label[0][1:].islower()
        ):
            label[0] = label[0][0] + '.'
        label = ' '.join(label)
        
        if showWeight:
            label = ('{} (%s)' % label_weight_formatting).format(
                label, 
                label_weight_factor * node.getWeight()
            )
            
        return label

    def export_json(
        self, 
        indent=0, 
        showWeight=False, 
        label_generator=None,
        label_weight_factor=100.0,
        label_weight_formatting='{:.2f}%'
    ):
        
        taxon = self.getTaxon()
        
        # generate label of node
        if label_generator is None:
            label_generator = TaxonTreeNode.generate_label
        label = label_generator(
            node=self,
            showWeight=showWeight,
            label_weight_factor=label_weight_factor,
            label_weight_formatting=label_weight_formatting
        )
        
        # export node in JSON format
        indent_spaces = ' ' * indent
        indent_spaces2 = ' ' * (indent - 4)
        export = '{\n'
        export += indent_spaces + '"^o":"Node", \n'
        export += indent_spaces + '"id":{},\n'.format(taxon.getId())
        export += indent_spaces + '"name":"{}",\n'.format(label)
        export += indent_spaces + '"data":{"weight":%d, "original_weight:%f"},\n' % (int(self.getWeight() * 255), self.getWeight())
        export += indent_spaces + '"children": ['
        export += ', '.join(
            child.export_json(
                indent=indent + 4, 
                showWeight=showWeight,
                label_weight_factor=label_weight_factor,
                label_weight_formatting=label_weight_formatting
            ) 
            for child in sorted(self.getChildren())
        )
        export += ']\n'
        export += indent_spaces2 + '}'
        
        # return JSON formatted node
        return export
    
    def export_table(
        self, 
        maxdepth=None, 
        minrank=None,
        showSubtotals=True, 
        prefix=None
    ):
        
        if maxdepth is not None and not maxdepth:
            return []
        
        # prefix is a list that is appended to each record; this shows the data
        # duplication that happens when hierarchical data structure is sliced 
        # and diced to convert it into a tabular data structure
        if prefix is None:
            prefix = []
            
        # append taxon to the already established prefix
        taxon = self.getTaxon()
        prefix += [
            taxon.getId(), 
            taxon.getName(), 
            taxon.getRank(), 
            self.getWeight()
        ]
        
        # recursively convert subhierarchies into tables
        if minrank != taxon.getRank():        
            recursive_records = [
                child.export_table(
                    maxdepth if maxdepth is None else maxdepth - 1,
                    minrank,
                    showSubtotals
                ) 
                for child in sorted(self.getChildren())
            ]
        else:
            recursive_records = []
        
        # flatten list of tables of records (2 levels) into list of records 
        output, child_weights = [], 0.0
        for table in recursive_records:
            for record in table:
                if maxdepth > 1:
                    output.append(prefix + record)
                    # exclude subtotals (rank -2), but include unidentified 
                    # (rank 0) and unspecified (rank -1) 
                    if record[-4] != -2:
                        child_weights += record[-1]

        # add weight of unspecified segment of node
        node_weight = self.getWeight()
        eps = 1e-6
        #if eps < childweights < node_weight - eps:
        if eps < child_weights:
            output.append(prefix + [-1, 'unspecified', 'no rank', node_weight - child_weights])
        elif not showSubtotals:
            output.append(prefix)
            
        # add node weight as subtotal if requested                
        if showSubtotals:
            if (
                prefix[0] >= 0
                and len(self.getChildren()) > 0
                and maxdepth > 1
                and minrank != taxon.getRank() 
            ):
                output.append(prefix + [-2, 'subtotal', 'no rank', node_weight])
            else:
                output.append(prefix)
            
        return output
    
    def normalizeNode(self, weight_tree=None):
        
        if weight_tree is None:
            # normalize weights across the entire tree, so that the sum of all
            # weights equals one
            weight_tree = self.getWeight()
            
        for node in self.nodes():
            node.setWeight(node.getWeight() / weight_tree)
            
    def normalizeChildren(self):
        
        node_weight = self.getWeight()
        for child in self.getChildren():
            child_weight = child.normalizeChildren()
            child.setWeight(child_weight / node_weight)
        return node_weight

    def __getitem__(self, key):
        
        for child in self.getChildren():
            if (
                (isinstance(key, int) and child.getTaxon().getId() == key)
                or
                (isinstance(key, str) and child.getTaxon().getName() == key)
            ):
                return child
            
        raise KeyError('no child found with key {}'.format(key))
    
    def countBranches(self):
        
        count, children = 0, 0
        for child in self.getChildren():
            children += 1
            count += child.countBranches()
        if children > 1:
            count += 1
        return count
        
class TaxonTree:
    
    """
    >>> import os
    >>> lineages = read_lineages(open(os.path.join('data', 'id_protein10.csv'), 'r'))
    >>> tree1 = TaxonTree()
    >>> for peptide, lineage in lineages: tree1.addLineage(lineage)
    >>> print(tree1)
    `-- cellular organism (1, root) (weight=2.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            `-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
                `-- Bacteroidetes (976, phylum) (weight=1.0)
                    `-- Bacteroidia (200643, class) (weight=0.0)
                        `-- Bacteroidales (171549, order) (weight=1.0)
                            `-- Bacteroidaceae (815, family) (weight=0.0)
                                `-- Bacteroides (816, genus) (weight=0.0)
                                    `-- Bacteroides fragilis (817, species) (weight=2.0)
                                    
    >>> tree1.redistributeWeightsTopDown()
    6.0
    >>> print(tree1)
    `-- cellular organism (1, root) (weight=0.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            `-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
                `-- Bacteroidetes (976, phylum) (weight=0.0)
                    `-- Bacteroidia (200643, class) (weight=0.0)
                        `-- Bacteroidales (171549, order) (weight=0.0)
                            `-- Bacteroidaceae (815, family) (weight=0.0)
                                `-- Bacteroides (816, genus) (weight=0.0)
                                    `-- Bacteroides fragilis (817, species) (weight=6.0)

    >>> tree1.normalizeWeights()
    >>> print(tree1)
    `-- cellular organism (1, root) (weight=0.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            `-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
                `-- Bacteroidetes (976, phylum) (weight=0.0)
                    `-- Bacteroidia (200643, class) (weight=0.0)
                        `-- Bacteroidales (171549, order) (weight=0.0)
                            `-- Bacteroidaceae (815, family) (weight=0.0)
                                `-- Bacteroides (816, genus) (weight=0.0)
                                    `-- Bacteroides fragilis (817, species) (weight=1.0)
                                    
    >>> tree1.percolateWeightsBottomUp()
    >>> print(tree1)
    `-- cellular organism (1, root) (weight=1.0)
        `-- Bacteria (2, superkingdom) (weight=1.0)
            `-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=1.0)
                `-- Bacteroidetes (976, phylum) (weight=1.0)
                    `-- Bacteroidia (200643, class) (weight=1.0)
                        `-- Bacteroidales (171549, order) (weight=1.0)
                            `-- Bacteroidaceae (815, family) (weight=1.0)
                                `-- Bacteroides (816, genus) (weight=1.0)
                                    `-- Bacteroides fragilis (817, species) (weight=1.0)

    >>> lineages = read_lineages(open(os.path.join('data', 'id_protein27.csv'), 'r'))
    >>> tree2 = TaxonTree()
    >>> for peptide, lineage in lineages: tree2.addLineage(lineage)
    >>> print(tree2)
    `-- cellular organism (1, root) (weight=3.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            |-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
            |   `-- Bacteroidetes (976, phylum) (weight=1.0)
            |       `-- Bacteroidia (200643, class) (weight=0.0)
            |           `-- Bacteroidales (171549, order) (weight=1.0)
            `-- Proteobacteria (1224, phylum) (weight=0.0)
                `-- Betaproteobacteria (28216, class) (weight=0.0)
                    `-- Burkholderiales (80840, order) (weight=0.0)
                        `-- Burkholderiaceae (119060, family) (weight=0.0)
                            `-- Burkholderia (32008, genus) (weight=0.0)
                                `-- Burkholderia cepacia complex (87882, species group) (weight=0.0)
                                    `-- Burkholderia cenocepacia (95486, species) (weight=1.0)

    >>> node = tree2.getNode(Taxon(32008, 'Burkholderia', 'genus'))
    >>> node.getTaxon()
    Taxon(32008, 'Burkholderia', 'genus')
    >>> for leaf in node.leafs(): print(repr(leaf))
    Burkholderia cenocepacia (95486, species) [weight=1.0]

    >>> for node in tree2.nodes(): print(repr(node))
    cellular organism (1, root) [weight=3.0]
    Bacteria (2, superkingdom) [weight=0.0]
    Bacteroidetes/Chlorobi group (68336, superphylum) [weight=0.0]
    Bacteroidetes (976, phylum) [weight=1.0]
    Bacteroidia (200643, class) [weight=0.0]
    Bacteroidales (171549, order) [weight=1.0]
    Proteobacteria (1224, phylum) [weight=0.0]
    Betaproteobacteria (28216, class) [weight=0.0]
    Burkholderiales (80840, order) [weight=0.0]
    Burkholderiaceae (119060, family) [weight=0.0]
    Burkholderia (32008, genus) [weight=0.0]
    Burkholderia cepacia complex (87882, species group) [weight=0.0]
    Burkholderia cenocepacia (95486, species) [weight=1.0]
    >>> for leaf in tree2.leafs(): print(repr(leaf))
    Bacteroidales (171549, order) [weight=1.0]
    Burkholderia cenocepacia (95486, species) [weight=1.0]
    
    >>> tree2.getTotalWeight()
    6.0
    >>> tree2.getTotalWeight(False)
    3.0
    
    >>> tree2.redistributeWeightsTopDown()
    6.0
    >>> print(tree2)
    `-- cellular organism (1, root) (weight=0.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            |-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
            |   `-- Bacteroidetes (976, phylum) (weight=0.0)
            |       `-- Bacteroidia (200643, class) (weight=0.0)
            |           `-- Bacteroidales (171549, order) (weight=4.0)
            `-- Proteobacteria (1224, phylum) (weight=0.0)
                `-- Betaproteobacteria (28216, class) (weight=0.0)
                    `-- Burkholderiales (80840, order) (weight=0.0)
                        `-- Burkholderiaceae (119060, family) (weight=0.0)
                            `-- Burkholderia (32008, genus) (weight=0.0)
                                `-- Burkholderia cepacia complex (87882, species group) (weight=0.0)
                                    `-- Burkholderia cenocepacia (95486, species) (weight=2.0)

    >>> tree2.normalizeWeights()
    >>> print(tree2)
    `-- cellular organism (1, root) (weight=0.0)
        `-- Bacteria (2, superkingdom) (weight=0.0)
            |-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.0)
            |   `-- Bacteroidetes (976, phylum) (weight=0.0)
            |       `-- Bacteroidia (200643, class) (weight=0.0)
            |           `-- Bacteroidales (171549, order) (weight=0.6666666666666666)
            `-- Proteobacteria (1224, phylum) (weight=0.0)
                `-- Betaproteobacteria (28216, class) (weight=0.0)
                    `-- Burkholderiales (80840, order) (weight=0.0)
                        `-- Burkholderiaceae (119060, family) (weight=0.0)
                            `-- Burkholderia (32008, genus) (weight=0.0)
                                `-- Burkholderia cepacia complex (87882, species group) (weight=0.0)
                                    `-- Burkholderia cenocepacia (95486, species) (weight=0.3333333333333333)

    >>> tree2.percolateWeightsBottomUp()
    >>> print(tree2)
    `-- cellular organism (1, root) (weight=1.0)
        `-- Bacteria (2, superkingdom) (weight=1.0)
            |-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.6666666666666666)
            |   `-- Bacteroidetes (976, phylum) (weight=0.6666666666666666)
            |       `-- Bacteroidia (200643, class) (weight=0.6666666666666666)
            |           `-- Bacteroidales (171549, order) (weight=0.6666666666666666)
            `-- Proteobacteria (1224, phylum) (weight=0.3333333333333333)
                `-- Betaproteobacteria (28216, class) (weight=0.3333333333333333)
                    `-- Burkholderiales (80840, order) (weight=0.3333333333333333)
                        `-- Burkholderiaceae (119060, family) (weight=0.3333333333333333)
                            `-- Burkholderia (32008, genus) (weight=0.3333333333333333)
                                `-- Burkholderia cepacia complex (87882, species group) (weight=0.3333333333333333)
                                    `-- Burkholderia cenocepacia (95486, species) (weight=0.3333333333333333)
                                    
    >>> merged = tree1 + tree2
    >>> print(merged)
    `-- cellular organism (1, root) (weight=1.0)
        `-- Bacteria (2, superkingdom) (weight=1.0)
            |-- Bacteroidetes/Chlorobi group (68336, superphylum) (weight=0.8333333333333333)
            |   `-- Bacteroidetes (976, phylum) (weight=0.8333333333333333)
            |       `-- Bacteroidia (200643, class) (weight=0.8333333333333333)
            |           `-- Bacteroidales (171549, order) (weight=0.8333333333333333)
            |               `-- Bacteroidaceae (815, family) (weight=0.5)
            |                   `-- Bacteroides (816, genus) (weight=0.5)
            |                       `-- Bacteroides fragilis (817, species) (weight=0.5)
            `-- Proteobacteria (1224, phylum) (weight=0.16666666666666666)
                `-- Betaproteobacteria (28216, class) (weight=0.16666666666666666)
                    `-- Burkholderiales (80840, order) (weight=0.16666666666666666)
                        `-- Burkholderiaceae (119060, family) (weight=0.16666666666666666)
                            `-- Burkholderia (32008, genus) (weight=0.16666666666666666)
                                `-- Burkholderia cepacia complex (87882, species group) (weight=0.16666666666666666)
                                    `-- Burkholderia cenocepacia (95486, species) (weight=0.16666666666666666)

    >>> merged['Bacteria']
    Bacteria (2, superkingdom) [weight=1.0]
    """
    
    def __init__(self):

        # initialize root of the tree        
        self.root = TaxonTreeNode(Taxon(1, 'cellular organism', 'root'))
        
        # initialize lookup table to lookup nodes in the tree
        self.id2node = {1:self.root}
        
    def getNode(self, taxon):
        
        # assure that node corresponding to given taxon occurs in the tree
        assert taxon.id in self.id2node, 'tree does not contain {}'.format(taxon)
        
        # lookup node in the tree
        return self.id2node[taxon.id]
        
    def addTaxon(self, taxon, parent_taxon=None):
        
        if taxon.id in self.id2node:
            
            # lookup taxon node if it already exists
            taxon_node = self.getNode(taxon)

        else:

            # lookup parent node            
            if parent_taxon is None:
                # set parent node to root if not parent given
                parent_taxon_node = self.root
                parent_taxon = parent_taxon_node.getTaxon()
            else:            
                # lookup node corresponding to parent taxon
                parent_taxon_node = self.getNode(parent_taxon)
            
            # add taxon node as child of parent node
            taxon_node = parent_taxon_node.addChild(taxon)
            
            # add taxon node to index
            self.id2node[taxon.id] = taxon_node
        
        return taxon_node
        
    def addLineage(self, lineage, weight=1.0):
        
        # lineage is always anchored at the root
        parent_taxon_node = self.root
        
        # link all taxa in the lineage to their parent taxon
        for taxon in lineage.getTaxa():
            if taxon is not None:
                parent_taxon_node = self.addTaxon(
                    taxon, 
                    parent_taxon_node.getTaxon()
                )
                
        # increment the weight of the child node
        parent_taxon_node.setWeight(parent_taxon_node.getWeight() + weight)
                
    def __str__(self):
        
        return str(self.root)
    
    def nodes(self, includeRoot=True):
        
        for node in self.root.nodes(includeRoot):
            yield node

    def leafs(self, includeRoot=True):
        
        for leaf in self.root.leafs(includeRoot):
            yield leaf
            
    def getTotalWeight(self, includeRoot=True):
        
        return self.root.getTotalWeight(includeRoot)
    
    def redistributeWeightsTopDown(self):
        
        # ignore weights at the root level, as they would be evenly distributed
        # over all leaf nodes anyway, which would have no effect after 
        # normalization
        # NOTE: now that I think about it, this is not entirely true; moreover
        #       it might change the weight of the tree itself
        # TODO: check the effect of removing the root weight
        #if self.root.getChildren():
        #    self.root.setWeight(0.0)
        
        return self.root.redistributeWeightsTopDown()

        # below is the initial implementation, which was top-down, but the
        # actual implementation (above) should be bottom-up (depth-first); I 
        # leave the actual implementation in, since this might be discussed in
        # the paper about this subject; after doing the implementation, it
        # seems as if both approaches always give the same result
        # TODO: it needs to be double-checked whether both implementations do
        #       produce the same result, and if so, which one is faster 
        
        # distribute non-zero weights of all internal nodes over the leaf nodes
        # in their subtree, according to the weights of the leaf nodes
        for node in self.nodes():
            if not node.isLeaf() and node.getWeight():
                
                # fetch node weight and then give node zero weight
                weight_node = node.getWeight()
                node.setWeight(0.0)
                
                # distribute node weight over its leaf nodes
                weight_subtree = node.getTotalWeight(
                    includeNode=False, 
                    leafsOnly=True
                )
                for leaf in node.leafs():
                    leaf.setWeight(
                        leaf.getWeight() * (1 + weight_node / weight_subtree)
                    )
                    
    def normalizeWeights(self, weight_tree=None):
        
        if weight_tree is None:
            # normalize weights across the entire tree, so that the sum of all
            # weights equals one
            weight_tree = self.root.getTotalWeight(
                includeNode=True, 
                leafsOnly=False
            )
            
        for node in self.nodes():
            node.setWeight(node.getWeight() / weight_tree)

    def percolateWeightsBottomUp(self):
        
        # compute weight of each internal node as the sum of the weights of its
        # child nodes; depth-first implementation for performance reasons
        self.root.percolateWeightsBottomUp()
        
    def __add__(self, other):
        
        return merge_trees({self, other})
    
    def export_json(
        self, 
        indent=0, 
        showWeight=False,
        label_weight_factor=100.0,
        label_weight_formatting='{:.2f}%'
    ):
        
        return self.root.export_json(
            indent=indent, 
            showWeight=showWeight,
            label_weight_factor=label_weight_factor,
            label_weight_formatting=label_weight_formatting
        )
        
    def export_csv(self, maxdepth=None, minrank=None):
        
        return '\n'.join(
            ','.join(
                str(field) 
                for field in record[4:]
            )
            for record in self.root.export_table(
                maxdepth + 1 if maxdepth is not None else None, 
                minrank
            )
        )
        
        # below is the old implementation, which did not show an unspecified
        # child for the root level; I'll leave the implementation in for a 
        # reference

        # discard missing proteins at root level        
        rootWeight = 0.0
        for child in self.root.getChildren():
            rootWeight += child.getWeight()
            
        for child in self.root.getChildren():
            child.normalizeNode(rootWeight)
                    
        return '\n'.join(
            '\n'.join(
                ','.join(
                    str(field) 
                    for field in record
                ) 
                for record in child.export_table(maxdepth, minrank)
            )
            for child in sorted(self.root.getChildren())
        )
    
    def __getitem__(self, key):
        return self.root.__getitem__(key)

    def normalizeChildren(self):
        return self.root.normalizeChildren()
    
    def countBranches(self):
        return self.root.countBranches()

def merge_trees(trees):
    
    # make a new tree which we will build into the merger of all trees
    merge = TaxonTree()
                    
    # merger of all trees has the union of all taxa in the individual trees
    # with the sum of the weights of the taxa in the individual trees
    for tree in trees:
        
        # add weights of all root elements
        merge.root.setWeight(merge.root.getWeight() + tree.root.getWeight())
        
        # add inner nodes and their weight
        for node in tree.nodes():
            parent_taxon = node.getTaxon()
            for child in node.getChildren():
                taxon = child.getTaxon()
                taxon_node = merge.addTaxon(taxon, parent_taxon)
                taxon_node.setWeight(
                    taxon_node.getWeight() + child.getWeight()
                )
    
    # normalize weights by the weight of the root
    merge.normalizeWeights(weight_tree=merge.root.getWeight())

    # return the merged tree         
    return merge

if __name__ == '__main__':
    import doctest
    doctest.testmod()