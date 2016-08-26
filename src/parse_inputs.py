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


def get_categ_samples_dict(mapping_info_list, factor):
    """Turn the groupfile in to a dictionary frob factor labels to sample ids
    """
    labels = mapping_info_list[0].strip().strip('#').split(run_config.DELIM)
    index_sampleid = labels.index('SampleID')
    index_categ = labels.index(factor)
    local_dict = dict()
    for l in mapping_info_list:
        if not l or l.strip()[0] == '#':
            continue
        key, val = map(l.strip().split(run_config.DELIM).__getitem__,
                       [index_categ, index_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict


def samples(mapping_dict, group):
    """Get the sample ids for the given group
"""
    return [otu for g in group for otu in mapping_dict[g]]


def num_present(vals, min_abundance):
    """Get the number of values in vals that are greater than min_abundance"""
    return sum([v > min_abundance for v in vals])
