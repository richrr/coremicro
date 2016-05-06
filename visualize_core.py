from ete3 import Tree, TreeStyle, NodeStyle, faces, TextFace
import re

# Feature name to mark nodes in the core OTUs
IN_CORE_FEATURE_NAME = 'in_core'
IN_OUT_FEATURE_NAME = 'in_out'
TAXA_FEATURE_NAME = 'taxa'

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
        if hasattr(node, IN_CORE_FEATURE_NAME):
            node.set_style(in_otu_style)
            faces.add_face_to_node(
                TextFace(getattr(node, TAXA_FEATURE_NAME)[-1], fgcolor='Red'),
                node, column=0)

        else:
            faces.add_face_to_node(
                TextFace(getattr(node, TAXA_FEATURE_NAME)[-1]),
                node, column=0)


# Style for exported trees
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_scale = False
ts.layout_fn = layout


def export_tree(tree, filename):
    '''
    exports the given tree to the given filename
    '''
    tree.render(filename, w=1080, tree_style=ts)


def mark_group(core_otus, tree, taxa_code_mapping, feature_name):
    '''
    Mark the given list of core otus in the given tree
    '''
    for otu in core_otus:
        mark_otu(otu, tree, taxa_code_mapping, feature_name)


def mark_otu(otu, tree, taxa_code_mapping, feature_name):
    '''
    Mark the given otu in the given tree
    '''
    otu_taxa = get_taxa(otu)
    code = taxa_code_mapping[otu_taxa]
    tree.search_nodes(name=code)[0].add_feature(feature_name, True)


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


def add_taxa(tree, code_taxa_mapping):
    '''
    add the taxa to the nodes
    '''
    for leaf in tree:
        leaf.add_feature(TAXA_FEATURE_NAME,
                         code_taxa_mapping[leaf.name.strip("'")])

if __name__ == '__main__':
    core = [
        # 100
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__',
        'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Alicyclobacillaceae;g__Alicyclobacillus;s__',
        'k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Cryomorphaceae;g__;s__',
        # 95
        'k__Bacteria;p__Bacteroidetes;c__[Saprospirae];o__[Saprospirales];f__Chitinophagaceae',
        'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
        # 90
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Phyllobacteriaceae;g__Mesorhizobium;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Nakamurellaceae;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Legionellales;f__;g__;s__'
        'k__Bacteria;p__;c__;o__;f__;g__;s__',
        'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__'
        # 85
        'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Nocardiaceae;g__Nocardia;s__',
        'k__Bacteria;p__Chloroflexi;c__TK10;o__;f__;g__;s__',
        'k__Bacteria;p__Chloroflexi;c__Anaerolineae;o__A31;f__S47;g__;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Intrasporangiaceae',
        'k__Bacteria;p__Elusimicrobia;c__Elusimicrobia;o__Elusimicrobiales;f__;g__;s__',
        'k__Bacteria;p__Chloroflexi;c__Chloroflexi;o__AKIW781;f__;g__;s__',
        'k__Bacteria;p__;c__;o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Novosphingobium;s__',
        'k__Bacteria;p__Chloroflexi;c__Ktedonobacteria;o__Ktedonobacterales;f__Ktedonobacteraceae;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Legionellales;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Phyllobacteriaceae;g__Mesorhizobium;s__',
        # 80
        'k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__[Weeksellaceae];g__Chryseobacterium;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Labrys;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
        'k__Bacteria;p__Bacteroidetes;c__[Saprospirae];o__[Saprospirales];f__Chitinophagaceae;g__Flavihumibacter;s__',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Methylophilales;f__Methylophilaceae;g__;s__',
        'k__Bacteria;p__Acidobacteria;c__BPC102;o__MVS-40;f__;g__;s__',
        'k__Bacteria;p__TM7;c__TM7-1;o__;f__;g__;s__',
        'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
        'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfuromonadales;f__Geobacteraceae;g__Geobacter;s__',
        'k__Bacteria;p__Acidobacteria;c__[Chloracidobacteria];o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__',
        'k__Bacteria;p__Chloroflexi;c__Anaerolineae;o__S0208;f__;g__;s__',
        'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Dyadobacter;s__',
        'k__Bacteria;p__Cyanobacteria;c__Chloroplast;o__Chlorophyta;f__Chlamydomonadaceae',
        'k__Bacteria;p__Cyanobacteria;c__Chloroplast;o__Chlorophyta;f__;g__;s__',
        'k__Bacteria;p__Chloroflexi;c__Chloroflexi;o__Herpetosiphonales;f__;g__;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Micromonosporaceae;g__Actinoplanes;s__',
        # 75
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Labrys;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Alteromonadales;f__OM60;g__;s__',
        'k__Bacteria;p__Acidobacteria;c__BPC102;o__MVS-40;f__;g__;s__',
        'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Pseudonocardiaceae',
        'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
        'k__Bacteria;p__Acidobacteria;c__BPC102;o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfuromonadales;f__Geobacteraceae;g__Geobacter;s__',
        'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Dyadobacter;s__',
        'k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__',
    ]
    out = [
        # 100
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
        'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
        # 95
        'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
        # 90
        'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
        'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
        # 85
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
        'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Myxococcales;f__Polyangiaceae;g__Sorangium;s__',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
        # 80
        'k__Bacteria;p__Armatimonadetes;c__[Fimbriimonadia];o__[Fimbriimonadales];f__[Fimbriimonadaceae];g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
        'k__Bacteria;p__Acidobacteria;c__Acidobacteria-6;o__;f__;g__;s__',
        'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
        'k__Bacteria;p__Acidobacteria;c__Acidobacteria-6;o__;f__;g__;s__',
    ]

    print 'importing'
    tree = load_gg_tree('97_otus.tree')
    print 'building map'
    mapping = build_taxonomy_code_map('97_otu_taxonomy.txt')
    print 'finding and marking core and out groups'
    mark_group(core, tree, mapping, IN_CORE_FEATURE_NAME)
    mark_group(out, tree, mapping, IN_OUT_FEATURE_NAME)
    core_nodes = tree.search_nodes(**{IN_CORE_FEATURE_NAME: True})
    out_nodes = tree.search_nodes(**{IN_OUT_FEATURE_NAME: True})
    print 'pruning...'
    tree.prune(core_nodes + out_nodes)

    add_taxa(tree, mapping)
    print tree.write(features=[])
    export_tree(tree, 'tree.png')
