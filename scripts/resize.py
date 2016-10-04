# Copyright 2016 Richard Rodrigues, Nyle Rodgers, Mark Williams, Virginia Tech
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
import argparse
from biom.parse import parse_biom_table
from biom.table import table_factory, SparseOTUTable
from numpy import array


def multiply_rows(table, multiply_by):
    """Returns a table with multiply_by times the number of rows of the given
    table. The new rows are identical to the rows of the original table, except
    that the ids are changed to make them unique. multiply_by should be an
    integer"""
    unused_data, sample_ids, sample_md = zip(*table.iterSamples())
    data, obs_ids, obs_md = zip(*table.iterObservations())

    data = [list(row) for row in data]
    data *= multiply_by
    obs_md *= multiply_by
    obs_ids = reduce(lambda a, b: a + b,
                     [['%s__%d' % (id, i) for id in obs_ids]
                      for i in range(multiply_by)])
    return table_factory(array(data), sample_ids, obs_ids, sample_md, obs_md,
                         constructor=SparseOTUTable)


def multiply_columns(table, multiply_by):
    """Returns a table with multiply_by times the number of columns of the given
    table. The new columns are identical to the columns of the original table,
    except that the ids are changed to make them unique. multiply_by should
    be an integer"""
    unused_data, sample_ids, sample_md = zip(*table.iterSamples())
    data, obs_ids, obs_md = zip(*table.iterObservations())

    data = [list(row) * multiply_by for row in data]
    sample_md *= multiply_by
    sample_ids = reduce(lambda a, b: a + b,
                        [['%s__%d' % (id, i) for id in sample_ids]
                         for i in range(multiply_by)])
    return table_factory(array(data), sample_ids, obs_ids, sample_md, obs_md,
                         constructor=SparseOTUTable)


def multiply_groupfile(groupfile, multiply_by):
    label_string = groupfile.split('\n')[0]
    labels = label_string.strip().strip('#').split('\t')
    id_index = labels.index('SampleID')
    rows = reduce(lambda a, b: a + b,
                  [[(lambda r: (r[0:id_index] + ['%s__%d' % (r[id_index], i)] +
                                r[id_index + 1:]))(row.split('\t'))
                    for row in groupfile.split('\n')[1:] if row.strip() != '']
                   for i in range(multiply_by)])
    return label_string + '\n' + '\n'.join(['\t'.join(row) for row in rows])


def setup_parser():
    parser = argparse.ArgumentParser(
        description='Enlarges the given data table by the given amount')
    parser.add_argument(
        'tablefile',
        help='The biom file containing the data table')
    parser.add_argument(
        'groupfile',
        help='The tsv file containing the grouping data')
    parser.add_argument(
        'table_outfile',
        help='the biom file to output the modified table to')
    parser.add_argument(
        'group_outfile',
        help='the tsv file to output the modified grouping info to')
    parser.add_argument(
        '-r', '--multiply_rows', type=int,  metavar='by', default=1,
        help='The number to multiply the number of rows by')
    parser.add_argument(
        '-c', '--multiply_columns', type=int, metavar='by', default=1,
        help='The number to multiply the number of columns by')
    return parser


if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    with open(args.tablefile) as f:
        table = parse_biom_table(f)
    with open(args.groupfile) as f:
        groupfile = f.read()
    table = multiply_rows(table, args.multiply_rows)
    table = multiply_columns(table, args.multiply_columns)
    groupfile = multiply_groupfile(groupfile, args.multiply_columns)
    with open(args.table_outfile, 'w') as f:
        f.write(table.getBiomFormatJsonString('COREMIC'))
    with open(args.group_outfile, 'w') as f:
        f.write(groupfile)
