# Copyright 2016, 2017 Richard Rodrigues, Nyle Rodgers, Mark Williams,
# Virginia Tech
#
# This file is part of Coremic.
#
# Coremic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Coremic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Coremic. If not, see <http://www.gnu.org/licenses/>.


from biom.parse import parse_biom_table
from biom.table import table_factory, SparseOTUTable

import run_config


def read_table(table_file):
    """Read in the input datafile and combine any doubled OTUs
    """
    parsed_table = parse_biom_table(table_file)
    sample_ids = parsed_table.SampleIds
    original_otus = len(parsed_table.ObservationIds)
    otu_data = dict()
    for vals, id, md in parsed_table.iterObservations():
        otu = md['taxonomy']
        if(isinstance(otu, type(list()))):
            otu = ';'.join(otu)
        if otu_data.get(otu) is not None:
            otu_data[otu] = vals + otu_data[otu]  # both are numpy arrays
        else:
            otu_data[otu] = vals

    observation_ids = otu_data.keys()
    data = [otu_data[o] for o in observation_ids]
    return table_factory(data, sample_ids, observation_ids,
                         constructor=SparseOTUTable), original_otus


def summarize_otu_data(table, mapping_dict, i_group, min_abundance):
    """Gives a summary of the data for each otu to use in processing
    """
    # column number of interest samples
    i_indexes = [i for i, id in enumerate(table.SampleIds)
                 if id in samples(mapping_dict, i_group)]
    interest = len(i_indexes)
    total = sum([len(mapping_dict[group]) for group in mapping_dict])
    otu_data = list()
    for vals, otu, md in table.iterObservations():
        i_present = len([v for i, v in enumerate(vals)
                         if i in i_indexes and v > min_abundance])
        present = len([v for v in vals if v > min_abundance])
        otu_data.append({
            'otu': otu,
            'total': total,
            'present': present,
            'absent': total - present,
            'interest': interest,
            'out': total - interest,
            'i_present': i_present,
            'i_absent': interest - i_present,
            'o_present': present - i_present,
            'o_absent': total - present - (interest - i_present),
        })
    return otu_data


def samples(mapping_dict, group):
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


def parse_inputs(params, mapping_file, data):
    """Validate that the given inputs are usable and parse them into usable
    formats"""
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
        filtered_data, original_otus = read_table(data)
        params['user_args'] += (
            'OTUs before combining duplicates: %s\n' +
            'OTUs after combining duplicates: %s\n\n\n') % (
                original_otus,
                len(filtered_data.ObservationIds)
            )
    except ValueError as e:
        errors_list.append('Datafile could not be read: %s' % e.message)

    if params['max_p'] < 0:
        errors_list.append('Maximum p-value can not be negative')
    elif params['max_p'] > 1:
        errors_list.append('Maximum p-value can not be greater than one')

    return (errors_list, mapping_dict, out_group, filtered_data)
