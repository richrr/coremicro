from biom.parse import parse_biom_table
from biom.table import table_factory, SparseOTUTable
from collections import namedtuple

import run_config


# Named Tuple to store information about the data in each row
OTU_Data = namedtuple('OTU_Data', ['n_interest', 'i_present',
                                   'total', 'present'])


def read_table(table_file):
    """Read in the input datafile and combine any doubled OTUs
    """
    parsed_table = parse_biom_table(table_file)
    sample_ids = parsed_table.SampleIds
    otu_data = dict()

    for vals, id, md in parsed_table.iterObservations():
        otu = md['taxonomy']
        if otu_data.get(otu) is not None:
            otu_data[otu] = vals + otu_data[otu]  # both are numpy arrays
        else:
            otu_data[otu] = vals

    observation_ids = otu_data.keys()
    data = [otu_data[o] for o in observation_ids]
    return table_factory(data, sample_ids, observation_ids,
                         constructor=SparseOTUTable)


def get_data_summary(table, mapping_dict, i_group, min_abundance):
    """Gives a summary of the data for each otu to use in processing
    """
    i_samples = samples(mapping_dict, i_group)
    n_interest = len(i_samples)
    total = sum([len(mapping_dict[group]) for group in mapping_dict])
    i_vals = {
        otu: vals for vals, otu, md in table.filterSamples(
            lambda values, id, md: id in i_samples
        ).iterObservations()
    }
    vals = {otu: vals for vals, otu, md in table.iterObservations()}
    return {otu: OTU_Data(n_interest,
                          num_present(i_vals[otu], min_abundance),
                          total,
                          num_present(vals[otu], min_abundance))
            for otu in i_vals.keys()}


def samples(mapping_dict, group):
    """Get the sample ids for the given group"""
    return [otu for g in group for otu in mapping_dict[g]]


def num_present(vals, min_abundance):
    """Get the number of values in vals that are greater than min_abundance"""
    return sum([v > min_abundance for v in vals])


def get_categ_samples_dict(mapping_info_list, factor):
    """Turn the groupfile in to a dictionary from factor labels to sample ids
    """
    # It is assumed that the first line is the header
    labels = mapping_info_list[0].strip().strip('#').split(run_config.DELIM)
    index_sampleid = labels.index('SampleID')
    index_categ = labels.index(factor)
    local_dict = dict()
    for l in mapping_info_list:
        if not l or l.strip()[0] == '#':  # Line is blank or comment
            continue
        key, val = map(l.strip().split(run_config.DELIM).__getitem__,
                       [index_categ, index_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict


def validate_inputs(params, inputs):
    """Validate that the given inputs are usable"""
    mapping_file = inputs['mapping_file']
    factor = params['factor']
    group = params['group']

    labels = mapping_file[0].strip().strip('#').split(run_config.DELIM)

    errors_list = list()

    if 'SampleID' not in labels:
        errors_list.append('"SampleID" not in the headers of the sample ' +
                           '<-> group info file')
    if factor not in labels:
        errors_list.append(('"%s" not in the headers of the sample <-> ' +
                            'group info file') % factor)
    mapping_dict = get_categ_samples_dict(mapping_file, factor)
    for l in group:
        if l not in mapping_dict.keys():
            errors_list.append(
                'Interest group label %s is not in groupfile' % l
            )

    try:
        out_group = mapping_dict.keys()
        [out_group.remove(l) for l in group]
    except ValueError:
        pass                    # already handled previously

    try:
        read_table(inputs['data'])
    except ValueError as e:
        errors_list.append('Datafile could not be read: %s' % e.message)

    if params['max_p'] < 0:
        errors_list.append('Maximum p-value can not be negative')
    elif params['max_p'] > 1:
        errors_list.append('Maximum p-value can not be greater than one')

    return (errors_list, mapping_dict, out_group)
