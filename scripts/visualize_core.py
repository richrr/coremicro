from ete3 import Tree, TreeStyle, faces, TextFace, RectFace
from colour import Color
import sys
import csv

# Feature name to mark nodes in the core OTUs
GROUP_FEATURE = 'group'
TAXA_FEATURE = 'taxa'
OTU_FEATURE = 'otu'
P_VAL_FEATURE = 'p_val'
CORRECTED_P_VAL_FEATURE = 'corrected_p_val'
THRESHOLD_FEATURE = 'threshold'

# The maximum corrected pvalue a significant otu might have
MAX_PVAL = 0.05                 # TODO: make this user changeable


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
            color = 'Black'
            if getattr(node, GROUP_FEATURE) == interest_group:
                pval = float(getattr(node, CORRECTED_P_VAL_FEATURE))
                pval_color = interpolate_color(signif_color, insignif_color,
                                               pval/MAX_PVAL)
                faces.add_face_to_node(RectFace(10, 10, fgcolor='Black',
                                                bgcolor=pval_color),
                                       node, column=1, aligned=True)
                faces.add_face_to_node(
                    TextFace('  {0:0.5f}  {1}'.format(
                        float(getattr(node, CORRECTED_P_VAL_FEATURE)),
                        getattr(node, OTU_FEATURE))),
                    node, column=2, aligned=True)
            else:  # node is in the out group
                color = 'Blue'
            faces.add_face_to_node(TextFace(getattr(node, TAXA_FEATURE)[-1],
                                            fgcolor=color),
                                   node, column=0)
    return layout


def interpolate_color(low_color, high_color, value):
    '''
    returns a color (hex code) for a color placed between low_color and
    high_color (both given as strings in a format that works with the colour
    library) scaled with value, which should be between zero and one; if
    value is zero, low_color will be returned, if value is one high_color
    will be returned. The scaling used is linear in HSL space.
    '''
    low = Color(low_color).hsl
    high = Color(high_color).hsl
    out = map((lambda h, l: (h - l) * value + l), high, low)
    out = Color(hue=out[0], saturation=out[1], luminance=out[2])
    return out.get_hex()


def export_tree(tree, filename, insignif_color, signif_color, i_group):
    '''
    exports the given tree to the given filename
    '''
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.layout_fn = generate_layout(insignif_color, signif_color, i_group)
    tree.render(filename, h=1080, tree_style=ts)


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


def add_group(otus, tree, group, threshold):
    '''
    Makes a tree from the given OTUs
    '''
    for otu in otus:
        if otu['threshold'] >= threshold:
            add_otu(otu, tree, group)
    return tree


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


if __name__ == '__main__':
    if len(sys.argv) not in [4, 5]:
        print("""Usage
visualize_core.py <interest_core.tsv> <out_core.tsv> <output.png> [threshold]

<interest_core.tsv> --- The file containing the core of the interest group
<out_core.tsv> --- The file containing the core of the out group
<output.png> --- The image file to output to
[threshold] --- (optional) The minimum threshold to include. 0 by default
""")
    i_core = parse_output(sys.argv[1])
    o_core = parse_output(sys.argv[2])
    if len(sys.argv) == 5:
        threshold = int(sys.argv[4])
    else:
        threshold = 0

    print 'building tree'
    tree = Tree()
    add_group(i_core, tree, 'i', threshold)
    add_group(o_core, tree, 'o', threshold)

    print 'exporting tree'
    c_signif = '#00441b'
    c_insignif = '#f7fcfd'
    export_tree(tree, sys.argv[3], c_signif, c_insignif, 'i')
