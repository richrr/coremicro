from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace

# Feature name to mark nodes in the core OTUs
IN_CORE_OTU_FEATURE_NAME = 'in_core'

# Style for nodes in the core OTUs
in_otu_style = NodeStyle()
in_otu_style['fgcolor'] = 'Red'
in_otu_style['hz_line_color'] = 'Red'


def layout(node):
    '''
    The layout function for exporting treees
    '''
    # Label the leaves
    if node.is_leaf():
        # Mark the core OTUs seprately
        if hasattr(node, IN_CORE_OTU_FEATURE_NAME):
            node.set_style(in_otu_style)
            faces.add_face_to_node(AttrFace('name', fgcolor='Red'), node,
                                   column=0)
        else:
            faces.add_face_to_node(AttrFace('name'), node, column=0)


# Style for exported trees
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_scale = False
ts.layout_fn = layout


def load_tree(filename, format=1):
    '''
    load and return a tree from the specified file
    '''
    return Tree(filename, format=format)


def export_tree(tree, filename):
    '''
    exports the given tree to the given filename
    '''
    tree.render(filename, w=1080, tree_style=ts)


def mark_core(core_otus, tree):
    '''
    Mark the given list of core otus in the given tree
    '''
    for otu in core_otus:
        mark_otu(otu, tree)


def mark_otu(otu, tree):
    '''
    Mark the given otu in the given tree
    '''
    otu_taxa = get_otu_taxa(otu)
    ptr = tree
    for taxon in otu_taxa:
        if not taxon:
            break
        child = find_immediate_child(ptr, taxon)
        if not child:           # The taxon isn't in the tree
            break
        ptr = child
    ptr.add_feature(IN_CORE_OTU_FEATURE_NAME, True)


def get_otu_taxa(name):
    '''
    Split up a name of format:
    k__[Kingdom][;p__[Phylum][;c__[Class][;o__[Order][;f__[Family][;g__[Genus][;s__[Species]]]]]]]

    return a list of the taxa of the otu
    '''
    return [c.strip()[3:] for c in name.split(';')]


def make_tree(core_otus):
    '''
    Makes a tree from the given OTUs
    '''
    tree = Tree()
    for otu in core_otus:
        add_otu(otu, tree)
    return tree


def add_otu(otu, tree):
    ptr = tree
    otu_taxa = get_otu_taxa(otu)
    for taxon in otu_taxa:
        if not taxon:
            break
        child = find_immediate_child(ptr, taxon)
        if not child:
            ptr = ptr.add_child(name=taxon)
        else:
            ptr = child


def find_immediate_child(node, name):
    '''
    looks for an immediate child of the given node with the given name. If
    none exists returns None, otherwise returns the first such child.
    '''
    for child in node.children:
        if child.name == name:
            return child
    return None


def remove_long_root(tree):
    '''
    If the root has only one child, uses that child as the new root. This is
    repeated until the root has multiple children
    '''
    while len(tree.children) == 1:
        tree = tree.children[0]
    return tree


def generate_tree_from_list(filename):
    '''
    Generates a tree from a list of OTUs in a file with the format
    <number>    k__[<Kingdom>][;p__[<Phylum>][;c__[<Class>][;o__[<Order>][;f__[<Family>][;g__[<Genus>][;s__[<Species>]]]]]]]
    '''
    file = open(filename, 'r')
    data = file.read()
    file.close()
    names = [line.split('\t')[1] for line in data.split('\n') if line]
    return make_tree(names)


if __name__ == '__main__':
    # tree = load_tree('big_tree.nw')
    core = [
        'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Alteromonadales;f__OM60;g__;s__',
        'k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__;s__ ',
        'k__Bacteria;p__Acidobacteria;c__BPC102;o__MVS-40;f__;g__;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Pseudonocardiaceae',
        'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
        'k__Bacteria;p__Acidobacteria;c__BPC102;o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfuromonadales;f__Geobacteraceae;g__Geobacter;s__',
        'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Dyadobacter;s__',
        'k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__',
        'k__Bacteria;p__WS3;c__PRR-12;o__Sediment-1;f__CV106;g__;s__',
        'k__Bacteria;p__Cyanobacteria;c__Chloroplast;o__Chlorophyta;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__',
    ]
    tree = make_tree(core)
    export_tree(tree, 'core.png')