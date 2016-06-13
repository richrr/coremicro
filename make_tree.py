from ete3.coretype.tree import Tree  # , TreeStyle, faces, TextFace, RectFace
import re
import logging

# from colour import Color

# Feature name to mark nodes in the core OTUs
IN_CORE_FEATURE = 'in_core'
IN_OUT_FEATURE = 'in_out'
TAXA_FEATURE = 'taxa'
OTU_FEATURE = 'otu'
FREQ_FEATURE = 'freq'
P_VAL_FEATURE = 'p_val'
CORRECTED_P_VAL_FEATURE = 'corrected_p_val'
THRESHOLD_FEATURE = 'threshold'

# The maximum corrected pvalue a significant otu might have
MAX_PVAL = 0.05


def annotate_group(group, tree, taxa_code_mapping, group_feature):
    '''
    Annotates the given group of otus in the given tree.
    Returns a list of the nodes annotated
    '''
    annotated = []
    for sig in group:
        for otu in group[sig]:
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
        FREQ_FEATURE: otu['freq'],
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


def make_tree(core, out):
    logging.info('importing')
    tree = load_gg_tree('tree_files/97_otus.tree')
    logging.info('building map')
    mapping = build_taxonomy_code_map('tree_files/97_otu_taxonomy.txt')
    logging.info('finding and marking core and out groups')
    core_nodes = annotate_group(core, tree, mapping, IN_CORE_FEATURE)
    out_nodes = annotate_group(out, tree, mapping, IN_OUT_FEATURE)
    logging.info('pruning')
    tree.prune(core_nodes + out_nodes)

    return tree.write(features=[], format=0)
