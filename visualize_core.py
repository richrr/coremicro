from ete3 import Tree, TreeStyle, faces, TextFace, RectFace
import re
from colour import Color
import sys
import csv

# Feature name to mark nodes in the core OTUs
IN_CORE_FEATURE = 'in_core'
IN_OUT_FEATURE = 'in_out'
TAXA_FEATURE = 'taxa'
OTU_FEATURE = 'otu'
P_VAL_FEATURE = 'p_val'
CORRECTED_P_VAL_FEATURE = 'corrected_p_val'
THRESHOLD_FEATURE = 'threshold'

# The maximum corrected pvalue a significant otu might have
MAX_PVAL = 0.05


def generate_layout(signif_color, insignif_color):
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
            if hasattr(node, IN_CORE_FEATURE):
                if hasattr(node, IN_OUT_FEATURE):
                    # This node is in the core for both the in and out group
                    # In theory this should not be happening
                    color = 'Red'
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


def export_tree(tree, filename, insignif_color, signif_color):
    '''
    exports the given tree to the given filename
    '''
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False
    ts.layout_fn = generate_layout(insignif_color, signif_color)
    tree.render(filename, h=1080, tree_style=ts)


def annotate_group(group, tree, taxa_code_mapping, group_feature, threshold):
    '''
    Annotates the given group of otus in the given tree.
    Returns a list of the nodes annotated
    '''
    annotated = []
    for otu in group:
        if otu['threshold'] >= threshold:
            annotated.append(annotate_otu(otu, tree, taxa_code_mapping,
                                          group_feature))
    return annotated


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


def load_gg_tree(filename):
    '''
    Loads an annotated tree from the GreenGenes database. These trees use a
    format that is a bit different from what ete3 wants.
    '''
    f = open(filename, 'r')
    data = f.read()
    f.close()
    # change things like (5:1.2) to ('5':1.2),
    # so that the data can work with
    data = re.sub('([(),])(\\d+\\.?\d*):(\\d+\\.?\\d*)', "\\1\'\\2\':\\3",
                  data)
    # replace : characters inside quotes with |
    data = re.sub("([(),]'[^(),]+?):([^().]+?')", '\\1|\\2', data)
    # replace ; characters not at the end with !
    data = re.sub(';(.)', '!\\1', data)
    return Tree(data, format=1)


def build_taxonomy_code_map(filename):
    '''
    build a dictonary maping taxonomies (stored as tuples) with their numeric
    codes (stored as strings with quotes)
    '''
    tax_code_map = dict()
    f = open(filename, 'r')
    for line in f:
        code, taxonomy_string = line.strip().split('\t')
        taxonomy = get_taxa(taxonomy_string)
        tax_code_map[code] = taxonomy
        tax_code_map[taxonomy] = "'" + code + "'"
    return tax_code_map


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


if __name__ == '__main__':
    if len(sys.argv) not in [4, 5]:
        print("""Usage
        visualize_core.py <core.tsv> <out.tsv> <out.png> [threshold]
""")
    core = parse_output(sys.argv[1])
    out = parse_output(sys.argv[2])
    if len(sys.argv) == 5:
        threshold = int(sys.argv[4])
    else:
        threshold = 0

    print 'importing'
    tree = load_gg_tree('tree_files/97_otus.tree')
    print 'building map'
    mapping = build_taxonomy_code_map('tree_files/97_otu_taxonomy.txt')
    print 'finding and marking core and out groups'
    core_nodes = annotate_group(core, tree, mapping, IN_CORE_FEATURE,
                                threshold)
    out_nodes = annotate_group(out, tree, mapping, IN_OUT_FEATURE,
                               threshold)

    print 'pruning...'
    tree.prune(core_nodes + out_nodes)

    # tree = Tree('pruned.nh')
    c_signif = '#00441b'
    c_insignif = '#f7fcfd'
    export_tree(tree, sys.argv[3], c_signif, c_insignif)
