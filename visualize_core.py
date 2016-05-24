from ete3 import Tree, TreeStyle, faces, TextFace, RectFace
import re
from colour import Color

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


def annotate_group(group, tree, taxa_code_mapping, group_feature):
    '''
    Annotates the given group of otus in the given tree.
    Returns a list of the nodes annotated
    '''
    annotated = []
    for otu in group:
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


if __name__ == '__main__':
    core = [
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Labrys;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': 0.04,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Alteromonadales;f__OM60;g__;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': 0.04,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__BPC102;o__MVS-40;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.024,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Pseudonocardiaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__BPC102;o__;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfuromonadales;f__Geobacteraceae;g__Geobacter;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Dyadobacter;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00342857142857,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.024,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.024,
            'threshold': 75
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__[Weeksellaceae];g__Chryseobacterium;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Labrys;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__[Saprospirae];o__[Saprospirales];f__Chitinophagaceae;g__Flavihumibacter;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': 0.04,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Methylophilales;f__Methylophilaceae;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__BPC102;o__MVS-40;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__TM7;c__TM7-1;o__;f__;g__;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': 0.04,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfuromonadales;f__Geobacteraceae;g__Geobacter;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__[Chloracidobacteria];o__;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Sphingobium;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__Anaerolineae;o__S0208;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Dyadobacter;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Cyanobacteria;c__Chloroplast;o__Chlorophyta;f__Chlamydomonadaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': 0.00485714285714,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Cyanobacteria;c__Chloroplast;o__Chlorophyta;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__Chloroflexi;o__Herpetosiphonales;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Micromonosporaceae;g__Actinoplanes;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': 0.0226666666667,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Nocardiaceae;g__Nocardia;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__TK10;o__;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__Anaerolineae;o__A31;f__S47;g__;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Intrasporangiaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Elusimicrobia;c__Elusimicrobia;o__Elusimicrobiales;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__Chloroflexi;o__AKIW781;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__;c__;o__;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .0266666666667,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Novosphingobium;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Chloroflexi;c__Ktedonobacteria;o__Ktedonobacterales;f__Ktedonobacteraceae;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00342857142857,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Legionellales;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .0266666666667,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Phyllobacteriaceae;g__Mesorhizobium;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00333333333333,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__;c__;o__;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00333333333333,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Legionellales;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .02,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Nakamurellaceae;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .02,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Phyllobacteriaceae;g__Mesorhizobium;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .00333333333333,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Planctomycetes;c__Planctomycetia;o__B97;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__[Saprospirae];o__[Saprospirales];f__Chitinophagaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Lysobacter;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .006,
            'threshold': 100,
        },
        {
            'otu': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Alicyclobacillaceae;g__Alicyclobacillus;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 100,
        },
        {
            'otu': 'k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Cryomorphaceae;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .03,
            'threshold': 100,
        },
    ]
    out = [
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__Acidobacteria-6;o__;f__;g__;s__',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .02,
            'threshold': 75,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Acidobacteria;c__Acidobacteria-6;o__;f__;g__;s__',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .008,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
            'freq': 1,
            'pval': 0.02,
            'corrected_pval': .04,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Armatimonadetes;c__[Fimbriimonadia];o__[Fimbriimonadales];f__[Fimbriimonadaceae];g__;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 80,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .003,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Myxococcales;f__Polyangiaceae;g__Sorangium;s__',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .003,
            'threshold': 85,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 90,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__A21b',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales;f__Solirubrobacteraceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .002,
            'threshold': 95,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Bradyrhizobiaceae;g__Bradyrhizobium',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .003,
            'threshold': 100,
        },
        {
            'otu': 'k__Bacteria;p__Actinobacteria;c__Thermoleophilia;o__Solirubrobacterales',
            'freq': 2,
            'pval': 0.04,
            'corrected_pval': .04,
            'threshold': 100,
        },
        {
            'otu': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Nitrosomonadales;f__Nitrosomonadaceae',
            'freq': 0,
            'pval': 0.0,
            'corrected_pval': .003,
            'threshold': 100,
        },

    ]

    print 'importing'
    tree = load_gg_tree('97_otus.tree')
    print 'building map'
    mapping = build_taxonomy_code_map('97_otu_taxonomy.txt')
    print 'finding and marking core and out groups'
    core_nodes = annotate_group(core, tree, mapping, IN_CORE_FEATURE)
    out_nodes = annotate_group(out, tree, mapping, IN_OUT_FEATURE)
    print 'pruning...'
    tree.prune(core_nodes + out_nodes)

    print tree.write(features=[])

    # tree = Tree('pruned.nh')
    c_signif = '#00441b'
    c_insignif = '#f7fcfd'
    export_tree(tree, 'tree.png', c_signif, c_insignif)
