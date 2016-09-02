from ete3 import Tree, TreeStyle, faces, TextFace, RectFace
from colour import Color
import sys
import csv
import argparse

# Feature name to mark nodes in the core OTUs
GROUP_FEATURE = 'group'
TAXA_FEATURE = 'taxa'
OTU_FEATURE = 'otu'
P_VAL_FEATURE = 'p_val'
CORRECTED_P_VAL_FEATURE = 'corrected_p_val'
THRESHOLD_FEATURE = 'threshold'


def generate_layout(signif_color, insignif_color, interest_group):
    '''
    Generates a layout function given a color for more insignigant otus
    and a color for more significant otus; a linear interpolation between
    these colors in HSB space will be used to generate intermediate colors
    from the p-value.
    '''
    def layout(node):
        '''
        The layout function for exporting treees
        '''
        # Label the leaves
        if node.is_leaf():
            faces.add_face_to_node(TextFace(
                getattr(node, TAXA_FEATURE)[-1],
                fgcolor=('blue'
                         if getattr(node, GROUP_FEATURE) == interest_group
                         else 'black')),
                                   node, column=0)
            pval = float(getattr(node, CORRECTED_P_VAL_FEATURE))
            pval_color = get_color(pval)
            faces.add_face_to_node(RectFace(10, 10, fgcolor='Black',
                                            bgcolor=pval_color),
                                   node, column=1, aligned=True)
            faces.add_face_to_node(
                TextFace('  {0:0.5f}  {1}'.format(
                    float(getattr(node, CORRECTED_P_VAL_FEATURE)),
                    getattr(node, OTU_FEATURE))),
                node, column=2, aligned=True)
    return layout


def get_color(value):
    if value < 0.001:
        return '#2ca25f'
    elif value < 0.01:
        return '#99d8c9'
    else:
        return '#e5f5f9'


def export_tree(tree, filename, insignif_color, signif_color, i_group, width):
    '''
    exports the given tree to the given filename
    '''
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.layout_fn = generate_layout(insignif_color, signif_color, i_group)
    tree.render(filename, w=width, tree_style=ts)


def annotate_otu(otu, tree, taxa_code_mapping, group_feature):
    '''
    Annotate the given otu in the given tree. Returns the node annotated
    '''
    otu_taxa = get_taxa(otu['otu'])
    code = taxa_code_mapping[otu_taxa]
    node = tree.search_nodes(name=code)[0]
    node.add_features(**{
        group_feature: True,
        TAXA_FEATURE: otu_taxa,
        OTU_FEATURE: otu['otu'],
        P_VAL_FEATURE: otu['pval'],
        CORRECTED_P_VAL_FEATURE: otu['corrected_pval'],
        THRESHOLD_FEATURE: otu['threshold'],
    })
    return node


def get_taxa(name):
    '''
    Split up a name of format:
    k__[Kingdom][;p__[Phylum][;c__[Class][;o__[Order][;f__[Family][;g__[Genus][;s__[Species]]]]]]]

    return a tuple of the taxa of the otu
    '''
    taxa = []
    for c in name.split(';'):
        c = c.strip()[3:]
        if c:
            taxa.append(c)
        else:
            break
    return tuple(taxa)


def parse_output(filename):
    core = []
    with open(filename) as tsv:
        reader = csv.reader(tsv, delimiter='\t')
        for row in reader:
            if row[0] == '' or row[0] == 'OTU':
                continue
            core.append({
                'otu': row[0],
                'pval': float(row[1]),
                'corrected_pval': float(row[2]),
                'threshold': int(row[3]),
            })
    return core


def add_group(otus, tree, group, min_threshold, max_threshold):
    '''
    Makes a tree from the given OTUs
    '''
    added_otus = 0
    for otu in otus:
        if otu['threshold'] >= min_threshold and\
           otu['threshold'] <= max_threshold:
            add_otu(otu, tree, group)
            added_otus += 1
    return added_otus


def find_immediate_child(node, name):
    '''
    looks for an immediate child of the given node with the given name. If
    none exists returns None, otherwise returns the first such child.
    '''
    for child in node.children:
        if child.name == name:
            return child
    return None


def add_otu(otu, tree, group):
    ptr = tree
    otu_taxa = get_taxa(otu['otu'])
    for taxon in otu_taxa:
        child = find_immediate_child(ptr, taxon)
        if not child:
            ptr = ptr.add_child(name=taxon)
        else:
            ptr = child
    ptr.add_features(**{
        GROUP_FEATURE: group,
        TAXA_FEATURE: otu_taxa,
        OTU_FEATURE: otu['otu'],
        P_VAL_FEATURE: otu['pval'],
        CORRECTED_P_VAL_FEATURE: otu['corrected_pval'],
        THRESHOLD_FEATURE: otu['threshold'],
    })


def trim_top(tree):
    ptr = tree
    while len(ptr.children) == 1:
        ptr = ptr.children[0]
    return ptr


def setup_parser():
    parser = argparse.ArgumentParser(
        description='Generate a tree visualization of the output from coremic')
    parser.add_argument(
        'interest_core',
        help='The tsv file containing the core of the interest group')
    parser.add_argument(
        'out_core',
        help='The tsv file containing the core of the out group')
    parser.add_argument(
        'output',
        help='The image file to output the result to')
    parser.add_argument(
        '-n', '--min_threshold', type=float, default=0,
        help='The minimum threshold to include in the tree. Defaults to zero.')
    parser.add_argument(
        '-x', '--max_threshold', type=float, default=100,
        help='The maximum threshold to include in the tree. Defaults to 100.')
    parser.add_argument(
        '-r', '--horizontal_resolution', type=int, default=1080,
        help='The horizontal resolution of the output image')
    return parser


if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    i_core = parse_output(args.interest_core)
    o_core = parse_output(args.out_core)
    tree = Tree()
    if not add_group(i_core, tree, 'i', args.min_threshold,
                     args.max_threshold)\
       and not add_group(o_core, tree, 'o', args.min_threshold,
                         args.max_threshold):
        print('No OTUs in specified threshold range!')
        sys.exit(1)
    tree = trim_top(tree)
    c_signif = '#00441b'
    c_insignif = '#f7fcfd'
    export_tree(tree, args.output, c_signif, c_insignif, 'i',
                args.horizontal_resolution)
