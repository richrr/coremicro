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
import argparse
import logging
from core.process_data import process, format_results
from core.parse_inputs import parse_inputs


def setup_parser():
    parser = argparse.ArgumentParser(
        description=('Runs Coremic analysis to find the core microbiome of ',
                     'the given interest group'))
    parser.add_argument('datafile', help=('BIOM file containing the data to ',
                                          'be analyzed'))
    parser.add_argument('groupfile', help=('Tab delimited file mapping each ',
                                           'SampleID to a group'))
    parser.add_argument('factor', help=('The factor with which the interest ',
                                        'group is identified. This should be ',
                                        'a column head in the groupfile'))
    parser.add_argument('group', help=('The value in the factor ',
                                       'column specifying the ',
                                       'interest group. Multiple ',
                                       'values may be specified ',
                                       'with commas separating them'))
    parser.add_argument('-p', '--max_p_val', type=float, metavar='pval',
                        default=0.05, help=('Maximum p-value to include in ',
                                            'the results; defaults to 0.05'))
    parser.add_argument('-t', '--min_presence', type=float, metavar='presence',
                        default=0.9, help=('Minimum fractional presence in ',
                                           'the interest group; defaults to ',
                                           '0.9 (90%)'))
    parser.add_argument('-o', '--max_out_presence', type=float,
                        metavar='presence',
                        default=1.0, help=('Maximum fractional presence in ',
                                           'the out group; defaults to ',
                                           '1.0 (100%)'))
    parser.add_argument('-a', '--max_absent_abundance', type=float,
                        metavar='abundance',
                        default=0.0, help=('Any abundance greater than this ',
                                           'value is considered to indicate ',
                                           'that the corresponding OTU is ',
                                           'present; defaults to 0'))
    parser.add_argument('-c', '--p_val_correction', default='b-h',
                        choices=['none', 'bf', 'bf-h', 'b-h'],
                        help=('The method to use for correcting for multiple ',
                              'testing; valid options are "none", "bf" ',
                              '(Bonferroni), "bf-h" (Bonferroni-Holm), and ',
                              '"b-h" (Benjamini-Hochberg, the default)'))
    parser.add_argument('-r', '--relative', dest='relative',
                        action='store_true',
                        help=('Convert the input datatable to relative ',
                              'abundance, rather than absolute abundance. ',
                              'I.E. abundance values are a percentage of ',
                              'each sample. If this is set then the maximum ',
                              'absent abundance should also be a relative ',
                              '(0 to 1) value.'))
    parser.set_defaults(relative=False)
    return parser


if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    with open(args.datafile) as f:
        datafile = f.read()
    with open(args.groupfile) as f:
        groupfile = f.read().split('\n')
    cfg = {
        'factor': args.factor,
        'group': map(lambda s: s.strip(), args.group.split(',')),
        'min_abundance': args.max_absent_abundance,
        'max_p': args.max_p_val,
        'min_frac': args.min_presence,
        'max_out_presence': args.max_out_presence,
        'p_val_adj': args.p_val_correction,
        'make_relative': args.relative,
    }
    errors_list, mapping_dict, out_group, filtered_data, original_otus \
        = parse_inputs(cfg, groupfile, datafile)
    if len(errors_list) > 0:
        logging.error(errors_list)
        exit(1)

    inputs = {
        'mapping_dict': mapping_dict,
        'filtered_data': filtered_data,
    }
    core = process(inputs, cfg)
    print(format_results(core, cfg))
    exit(0)
