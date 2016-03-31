#!/usr/bin/env python
# File created on 12 Jun 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from os.path import join
#from matplotlib import use
#use('Agg', warn=False)
#from pylab import xlim, ylim, xlabel, ylabel, plot, savefig
from numpy import linspace
#from cogent.util.misc import create_dir

from string import strip

from biom.parse import parse_biom_table
from biom.exception import TableException
#from qiime.util import parse_command_line_parameters, make_option
#from qiime.format import format_biom_table
from qiime.core_microbiome import filter_table_to_core
#from qiime.filter import sample_ids_from_metadata_description

#### check why qiime.filter would not access the methods but biom.parse can
### maybe because those qiime scripts needs a lot other scripts

'''
script_info = {}
script_info['brief_description'] = "Identify the core microbiome."
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("","Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in at least 50% of the samples. Write the list of core OTUs to a text file, and a new BIOM file containing only the core OTUs.","%prog -i otu_table.biom -o otu_table_core"))

script_info['script_usage'].append(("","Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in all of the samples in the 'Fast' treatment (as specified in the mapping file). Write the list of core OTUs to a text file.","%prog -i otu_table.biom -o otu_table_core_fast --mapping_fp map.txt --valid_states \"Treatment:Fast\""))


script_info['output_description']= ""
script_info['required_options'] = [\
 make_option('-i','--input_fp',type="existing_filepath",help='the input otu table in BIOM format'),
 make_option('-o','--output_dir',type="new_dirpath",help='directory to store output data'),
]
script_info['optional_options'] = [
 make_option('--max_fraction_for_core',type=float,
             help='the max fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: %default]',default=1.0),
 make_option('--min_fraction_for_core',type=float,
             help='the min fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: %default]',default=0.5),
 make_option('--num_fraction_for_core_steps',type=int,
             help='the number of evenly sizes steps to take between min_fraction_for_core and max_fraction_for_core [default: %default]',default=11),
 make_option('--otu_md',default='taxonomy',
             help='the otu metadata category to write to the output file [defualt: %default]'),\
 make_option('--mapping_fp',type='existing_filepath',
  help='mapping file path (for use with --valid_states) [default: %default]'),
 make_option('--valid_states',
  help='description of sample ids to retain (for use with --mapping_fp) [default: %default]')
]
script_info['version'] = __version__
'''


def format_biom_table(biom_table):
    """ Given a biom-format Table object, returns that Table as a BIOM string"""
    generated_by_str = "QIIME 1.8.0"
    return biom_table.getBiomFormatJsonString(generated_by_str)


def parse_mapping_file(lines, strip_quotes=True, suppress_stripping=False):
    """Parser for map file that relates samples to metadata.
    
    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    Result: list of lists of fields, incl. headers.
    """
    if hasattr(lines,"upper"):
        # Try opening if a string was passed
        try:
            lines = open(lines,'U')
        except IOError:
            raise QiimeParseError,\
             ("A string was passed that doesn't refer "
              "to an accessible filepath.")
        
    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"','')
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"','').strip()
    else:
        if suppress_stripping:
            # don't remove quotes or spaces
            strip_f = lambda x: x
        else:
            # remove spaces but not quotes
            strip_f = lambda x: x.strip()
    
    # Create lists to store the results
    mapping_data = []
    header = []
    comments = []
    
    # Begin iterating over lines
    for line in lines:
        line = strip_f(line)
        if not line or (suppress_stripping and not line.strip()):
            # skip blank lines when not stripping lines
            continue
        
        if line.startswith('#'):
            line = line[1:]
            if not header:
                header = line.strip().split('\t')
            else:
                comments.append(line)
        else:
            # Will add empty string to empty fields
            tmp_line = map(strip_f, line.split('\t'))
            if len(tmp_line)<len(header):
                tmp_line.extend(['']*(len(header)-len(tmp_line)))
            mapping_data.append(tmp_line)
    if not header:
        raise QiimeParseError, "No header line was found in mapping file."
    if not mapping_data:
        raise QiimeParseError, "No data found in mapping file."
    
    return mapping_data, header, comments

def parse_metadata_state_descriptions(state_string):
    """From string in format 'col1:good1,good2;col2:good1' return dict."""
    result = {}
    state_string = state_string.strip()
    if state_string:
        cols = map(strip, state_string.split(';'))
        for c in cols:
            # split on the first colon to account for category names with colons
            colname, vals = map(strip, c.split(':', 1))
            vals = map(strip, vals.split(','))
            result[colname] = set(vals)
    return result

def sample_ids_from_metadata_description(mapping_f,valid_states_str):
    """ Given a description of metadata, return the corresponding sample ids
    """
    map_data, map_header, map_comments = parse_mapping_file(mapping_f)
    valid_states = parse_metadata_state_descriptions(valid_states_str)
    sample_ids = get_sample_ids(map_data, map_header, valid_states)

    if len(sample_ids)<1:
        raise ValueError,"All samples have been filtered out for the criteria"+\
            " described in the valid states"

    return sample_ids

def get_sample_ids(map_data, map_header, states):
    """Takes col states in {col:[vals]} format.

    If val starts with !, exclude rather than include.
    
    Combines cols with and, states with or.

    For example, Study:Dog,Hand will return rows where Study is Dog or Hand;
    Study:Dog,Hand;BodySite:Palm,Stool will return rows where Study is Dog
    or Hand _and_ BodySite is Palm or Stool; Study:*,!Dog;BodySite:*,!Stool
    will return all rows except the ones where the Study is Dog or the BodySite
    is Stool.
    """
    
    name_to_col = dict([(s,map_header.index(s)) for s in states])
    good_ids = []
    for row in map_data:    #remember to exclude header
        include = True
        for s, vals in states.items():
            curr_state = row[name_to_col[s]]
            include = include and (curr_state in vals or '*' in vals) \
                and not '!'+curr_state in vals
        if include:        
            good_ids.append(row[0])
    return good_ids


def exec_core_microb_cmd(infile, output_dir, mapping_info_list, categ, group):
    #cmd = 'compute_core_microbiome.py -i %s -o %s --mapping_fp %s --valid_states "%s:%s"' %(infile, out_dir, new_mapping_file, categ, group)
    
    min_fraction_for_core = 0.9 # change this to 0.5 or 0.75 if needed
    max_fraction_for_core = 1.0
    num_fraction_for_core_steps = 3  # change this to 11 if it is 0.5 above; 6 if 0.75 above
    
    fractions_for_core = linspace(min_fraction_for_core,
                                  max_fraction_for_core,
                                  num_fraction_for_core_steps)

    otu_md = 'taxonomy'
    valid_states = "%s:%s" %(categ, group)
    
    #create_dir(output_dir)  ### this may not work

    if valid_states and mapping_info_list:
        sample_ids = sample_ids_from_metadata_description(mapping_info_list, valid_states)
        if len(sample_ids) < 1:
            return "--valid_states pattern didn't match any entries in mapping file: \"%s\"" % valid_states
    else:
        # get core across all samples if user doesn't specify a subset of the 
        # samples to work with
        sample_ids = None
    
    input_table = parse_biom_table(infile)
    
    otu_counts = []
    summary_figure_fp = join(output_dir,'core_otu_size.pdf')
    
    frac_thresh_core_OTUs_biom = dict()
    # for every fraction threshold, we get the:
    #   list of core OTUs (list)
    #   core biom table (string)
        
    for fraction_for_core in fractions_for_core:
        # build a string representation of the fraction as that gets used 
        # several times
        fraction_for_core_str = "%1.0f" % (fraction_for_core * 100.)
        
        # prep output files
        #output_fp = join(output_dir,'core_otus_%s.txt' % fraction_for_core_str)
        #output_table_fp = join(output_dir,'core_table_%s.biom' % fraction_for_core_str)
        output_f = list()  # this contains the list of core OTUs to be written eventually
    
        try:
            core_table = filter_table_to_core(input_table,
                                              sample_ids,
                                              fraction_for_core)
        except TableException:
            return "# No OTUs present in %s %% of samples." % fraction_for_core_str
            otu_counts.append(0)
            continue
    
        # write some header information to file
        if sample_ids == None:
            output_f.append("# Core OTUs across %s %% of samples.\n" % fraction_for_core_str)
        else:
            output_f.append(\
             "# Core OTUs across %s %% of samples matching the sample metadata pattern \"%s\":\n# %s\n" %\
              (fraction_for_core_str, valid_states,' '.join(sample_ids)))
    
        # write the otu id and corresponding metadata for all core otus
        otu_count = 0
        for value, id_, md in core_table.iterObservations():
            output_f.append('%s\t%s\n' % (id_,md[otu_md]))
            otu_count += 1

        # write the core biom table
        output_table_f = format_biom_table(core_table)
        
        frac_thresh_core_OTUs_biom[fraction_for_core_str] = [ output_f, output_table_f]
        
        # append the otu count to the list of counts
        otu_counts.append(otu_count)
    
    return {'fractions_for_core' : fractions_for_core, 'otu_counts' : otu_counts,\
         'frac_thresh_core_OTUs_biom' : frac_thresh_core_OTUs_biom }

