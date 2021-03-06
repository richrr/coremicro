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
from pval import getpval, correct_pvalues
from otu import Otu


def process(inputs, cfg):
    """Finds the core OTUs"""
    interest_ids = [otu for g in cfg['group']
                    for otu in inputs['mapping_dict'][g]]
    i_indexes = [i for i, id in enumerate(inputs['filtered_data'].SampleIds)
                 if id in interest_ids]
    otus = [Otu(vals, name, cfg['min_abundance'], i_indexes)
            for vals, name, md in inputs['filtered_data'].iterObservations()]
    pvals = list()
    for otu in otus:
        pval = getpval(otu)
        otu.pval = pval
        pvals.append(pval)
    for otu, corrected_pval in zip(otus,
                                   correct_pvalues(pvals, cfg['p_val_adj'])):
        otu.corrected_pval = corrected_pval
    # Filter down to the core
    return [otu for otu in otus
            if (otu.corrected_pval <= cfg['max_p'] and
                otu.interest_frac >= cfg['min_frac'] and
                otu.out_frac <= cfg['max_out_presence'])]


def format_results(res, cfg):
    """Format the result data as a tsv"""
    # Summary of the inputs given
    inputs = ['#Factor: ' + cfg['factor'],
              'Group: ' + ', '.join(cfg['group']),
              'Max Corrected p-val: %f' % cfg['max_p'],
              'Min Presence: %f' % cfg['min_frac'],
              'Max Out Presence: %f' % cfg['max_out_presence'],
              'Min Abundance: %f' % cfg['min_abundance'],
              'Correction Type: ' + cfg['p_val_adj']]
    # Header for results
    header = ['#OTU', 'Pval', 'Corrected Pval', 'Interest Group Presence',
              'Out Group Presence']
    # Combine inputs, header, and information from core OTUs
    return to_tsv([inputs, header] +
                  [[otu.name, otu.pval, otu.corrected_pval, otu.interest_frac,
                    otu.out_frac] for otu in list(sorted(res))])


def to_tsv(values):
    """Formats the given list of lists as a tsv string. str() is called on
    all items to convert them to strings"""
    return '\n'.join(map(lambda r: '\t'.join(map(str, r)), values))
