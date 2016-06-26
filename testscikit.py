
from biom.table import SparseOTUTable, DenseOTUTable, table_factory
import json
from numpy import zeros, asarray, uint32, float64
from string import strip


######################################################################
# READ FILE, ARG: FILENAME
# RETURN LIST
#####################################################################
def read_file(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    #print 'Read %d lines from file %s ' %(len(lines), filename)
    return lines

# the functions you need from biom.parse
def parse_biom_table(json_fh,constructor=None):
    """parses a biom format otu table into a rich otu table object

    input is an open filehandle or compatable object (e.g. list of lines)

    sparse/dense will be determined by "matrix_type" in biom file, and 
    either a SparseOTUTable or DenseOTUTable object will be returned
    note that sparse here refers to the compressed format of [row,col,count]
    dense refers to the full / standard matrix representations

    If try_light_parse is True, the light_parse_biom_sparse call will be 
    attempted. If that parse fails, the code will fall back to the regular
    BIOM parser.
    """
    table_str = ''.join(json_fh)

    t = parse_biom_table_str(table_str, constructor=constructor)
    return t

def parse_biom_table_str(json_str,constructor=None, data_pump=None):
    """Parses a JSON string of the Biom table into a rich table object.
   
    If constructor is none, the constructor is determined based on BIOM
    information

    data_pump is to allow the injection of a pre-parsed data object
    """
    json_table = json.loads(json_str)

    if constructor is None:
        f = BIOM_TYPES.get(json_table['type'].lower(), None)
    else:
        f = BIOM_TYPES.get(constructor._biom_type.lower(), None)

        # convert matrix data if the biom type doesn't match matrix type
        # of the table objects
        if constructor._biom_matrix_type != json_table['matrix_type'].lower():
            if json_table['matrix_type'] == 'dense':
                # dense -> sparse
                conv_data = []
                for row_idx,row in enumerate(json_table['data']):
                    for col_idx, value in enumerate(row):
                        if value == 0:
                            continue
                        conv_data.append([row_idx,col_idx,value])
                json_table['data'] = conv_data

            elif json_table['matrix_type'] == 'sparse':
                # sparse -> dense
                conv_data = zeros(json_table['shape'],dtype=float)
                for r,c,v in json_table['data']:
                    conv_data[r,c] = v
                json_table['data'] = [list(row) for row in conv_data]

            else:
                raise BiomParseException, "Unknown matrix_type"

    if f is None:
        raise BiomParseException, 'Unknown table type'

    return f(json_table, constructor, data_pump)

def convert_biom_to_table(biom_f, header_key=None, header_value=None, \
        md_format=None):
    """Convert a biom table to a contigency table"""
    table = parse_biom_table(biom_f)

    if md_format is None:
        md_format = biom_meta_to_string

    if table.ObservationMetadata is None:
        return table.delimitedSelf()
    
    if header_key in table.ObservationMetadata[0]:
        return table.delimitedSelf(header_key=header_key, 
                                       header_value=header_value,
                                       metadata_formatter=md_format)
    else:
        return table.delimitedSelf()


def pick_constructor(mat_type, table_type, constructor, valid_constructors):
    """Make sure constructor is sane, attempt to pick one if not specified

    Excepts valid_constructors to be a list in the order of 
    [SparseTable, DenseTable] in which the objects present must subclass the
    objects respectively (eg [SparseOTUTable, DenseOTUTable])

    We do not require the matrix type to be the same as the constructor if the 
    passed in constructor is not None. The motivation is that there are use
    cases for taking a table stored as dense but loaded as sparse.

    Will raise BiomParseError if input_mat_type appears wrong or if the 
    specified constructor appears to be incorrect
    """
    if constructor is None:
        if mat_type.lower() == 'sparse':
            constructor = valid_constructors[0]
        elif mat_type.lower() == 'dense':
            constructor = valid_constructors[1]
        else:
            raise BiomParseException, "Unknown matrix_type"

    if constructor._biom_type.lower() != table_type.lower():
        raise BiomParseException, "constructor must be a biom %s" % table_type

    return constructor


def parse_biom_otu_table(json_table, constructor=None, data_pump=None):
    """Parse a biom otu table type

    Constructor must have a _biom_type of "otu table"
    """
    table_type = 'otu table'
    mat_type = json_table['matrix_type']
    constructors = [SparseOTUTable, DenseOTUTable]
    constructor = pick_constructor(mat_type,table_type,constructor,constructors)

    sample_ids = [col['id'] for col in json_table['columns']]
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    obs_ids = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]
    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]

    if data_pump is None:
        table_obj = table_factory(json_table['data'], sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)
    else:
        table_obj = table_factory(data_pump, sample_ids, obs_ids, 
                                  sample_metadata, obs_metadata, 
                                  constructor=constructor, 
                                  shape=json_table['shape'], 
                                  dtype=dtype)

    return table_obj


# convert biom to table
def main_converter(input_fp, header_key, biom_to_classic_table=True):
    
    if input_fp is None:
        return ('Must specify an input file!')

    input_f = input_fp    
    output_metadata_id = header_key

    if biom_to_classic_table:
        try:
            return(convert_biom_to_table(input_f, header_key, output_metadata_id, lambda x: x)) # to avoid every character in taxonomy being split and joined by '; '
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"


# map table types -> parsing methods
BIOM_TYPES = {'otu table':parse_biom_otu_table}

MATRIX_ELEMENT_TYPE = {'int':int,'float':float,'unicode':unicode,
                      u'int':int,u'float':float,u'unicode':unicode}

QUOTE = '"'
JSON_OPEN = set(["[", "{"])
JSON_CLOSE = set(["]", "}"])
JSON_SKIP = set([" ","\t","\n",","])
JSON_START = set(["0","1","2","3","4","5","6","7","8","9","{","[",'"'])



# read a biom file as a list of lines
biom_list_of_lines = read_file("static/sample_data/otu_table_16s.biom")

# convert to classical table
otu_table = main_converter(biom_list_of_lines, header_key="taxonomy")
#print otu_table

# convert table to scikit input
