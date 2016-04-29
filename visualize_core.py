from ete3 import Tree, TreeStyle, NodeStyle

# The feature name to use to annotate if something is in the core
IN_CORE_OTU_FEATURE_NAME = 'in_core'


def load_tree(filename):
    '''
    load and return a tree from the specified file
    '''
    return Tree(filename, format=8)


def export_marked_tree(tree, filename):
    '''
    exports the given tree to the given filename, hilighting the marked
    core otus
    '''
    ts = TreeStyle()
    ts.show_leaf_name = True

    in_otu_style = NodeStyle()
    in_otu_style["vt_line_color"] = 'Red'
    in_otu_style["vt_line_width"] = 4
    in_otu_style["hz_line_color"] = 'Red'
    in_otu_style["hz_line_width"] = 4

    marked = tree.search_nodes(IN_CORE_OTU_FEATURE_NAME=True)
    for node in marked:
        node.set_Style(in_otu_style)

    tree.render(filename, w=400, tree_style=ts)  #


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
    for taxa in otu_taxa:
        if not taxa:
            break
        for child in tree.children:
            if child.name == taxa:
                tree = child
                break
    tree.add_feature(IN_CORE_OTU_FEATURE_NAME, True)


def get_otu_taxa(name):
    '''
    Split up a name of format:
    k__[Kingdom][;p__[Phylum][;c__[Class][;o__[Order][;f__[Family][;g__[Genus][;s__[Species]]]]]]]

    return a list of the taxa of the otu
    '''
    return [c[3:] for c in name.split(';')]


if __name__ == '__main__':
    tree = load_tree('tree.nh')
    # 'k__;p__;c__;o__;f__;g__;s__',
    mark_core([
        'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Alteromonadales;f__OM60;g__;s__'
    ], tree)
    export_marked_tree(tree, 'testoutput.png')
