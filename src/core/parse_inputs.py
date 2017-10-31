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
from numpy import array, append
from itertools import groupby

# groupfile delimitere
DELIM = '\t'


def read_table(table_file, normalize=False, do_quantile_normalize=False):
    """Read in the input datafile and combine any doubled OTUs
    """
    parsed_table = parse_biom_table(table_file)
    sample_ids = parsed_table.SampleIds
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
    if normalize:
        data = normalize_columns(data)
    if do_quantile_normalize:
        data = quantile_normalize(data)
    return table_factory(data, sample_ids, observation_ids,
                         constructor=SparseOTUTable)


def combine_tables(tables):
    """Combines multiple biom tables into a signle table, discarding any
    non-shared OTUs.
    """
    samples = [sample for sample_list in [table.SampleIds for table in tables]
               for sample in sample_list]
    duplicate_sample_indices = [
        index for indices in [[index for index, value in indices][1:]
                              for sample, indices
                              in groupby(sorted(enumerate(samples),
                                                key=lambda x: x[1]),
                                         lambda x: x[1])]
        for index in indices]
    otu_data = dict()
    for table in tables:
        for vals, otu, md in table.iterObservations():
            if otu_data.get(otu) is not None:
                otu_data[otu] = append(otu_data[otu], vals)
            else:
                otu_data[otu] = vals
    otus = [otu for otu in otu_data if len(otu_data[otu]) == len(samples)]
    if not otus:
        raise ValueError('No shared OTUs')
    data = [array([v for i, v in enumerate(otu_data[otu])
                   if i not in duplicate_sample_indices]) for otu in otus]
    samples = [v for i, v in enumerate(samples)
               if i not in duplicate_sample_indices]

    return table_factory(data, samples, otus, constructor=SparseOTUTable)


def normalize_columns(data):
    """Adjust the given two dimensional array so that the sum of the values in
    each column in one"""
    for i in range(len(data[0])):
        total = sum([row[i] for row in data])
        for row in data:
            row[i] = row[i] / float(total)
    return data


def quantile_normalize(data):
    """Quantile normalize the data"""
    columns = [[row[i] for row in data] for i in range(len(data[0]))]
    sorted_columns = [sorted(column) for column in columns]
    distribution = [(sum([column[i] for column in sorted_columns]) /
                     len(sorted_columns))
                    for i in range(len(sorted_columns[0]))]
    result = [
        [newval for newval, column_i
         in sorted([(distribution[dist_i], column_i)
                    for dist_i, (column_i, val)
                    in enumerate(sorted(enumerate(column),
                                        key=lambda x: x[1]))],
                   key=lambda x: x[1])]
        for column in columns]
    return [array([column[i] for column in result])
            for i in range(len(result[0]))]


def parse_groupfile(groupfile, factor):
    """Turn the groupfile in to a dictionary from factor labels to sample ids
    """
    # It is assumed that the first line is the header
    labels = groupfile[0].strip().strip('#').split(DELIM)
    if 'SampleID' not in labels:
        raise ValueError('"SampleID" not in the headers of the groupfile')
    if factor not in labels:
        raise ValueError('"%s" not in the headers of the groupfile' % factor)

    index_sampleid = labels.index('SampleID')
    index_categ = labels.index(factor)
    local_dict = dict()
    for l in groupfile:
        if not l or l.strip()[0] == '#':  # Line is blank or comment
            continue
        key, val = map(l.strip().split(DELIM).__getitem__,
                       [index_categ, index_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict


def parse_inputs(params, groupfile, datafiles):
    """Validate that the given inputs are usable and parse them into usable
    formats"""

    errors_list = list()

    try:
        mapping_dict = parse_groupfile(groupfile, params['factor'])
        for l in params['group']:
            if l not in mapping_dict.keys():
                errors_list.append(
                    'Interest group label %s is not in groupfile' % l)

        out_group = mapping_dict.keys()
        [out_group.remove(l) for l in params['group']]
    except ValueError as e:
        errors_list.append(e.message)
        mapping_dict = dict()

    try:
        filtered_data = combine_tables(map(
            lambda table: read_table(table, params['make_relative'],
                                     params['quantile_normalize']),
            datafiles))
    except ValueError as e:
        errors_list.append('Datafile could not be read: %s' % e.message)
        filtered_data = None

    if params['max_p'] < 0 or params['max_p'] > 1:
        errors_list.append('Maximum p-value must be in the range zero to one')
    if params['min_frac'] < 0 or params['min_frac'] > 1:
        errors_list.append('Minimum interest group presence must be in range' +
                           ' zero to one')
    if params['max_out_presence'] < 0 or params['max_out_presence'] > 1:
        errors_list.append('Maximum out-group presence must be in range' +
                           ' zero to one')
    # Should verify p-value-adjust is valid value

    return (errors_list, mapping_dict, out_group, filtered_data)
