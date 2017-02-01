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
from correct_p_values import correct_pvalues
from otu import Otu


def process(inputs, cfg):
    """Finds the core OTUs"""
    interest_ids = [otu for g in cfg['group']
                    for otu in inputs['mapping_dict'][g]]
    i_indexes = [i for i, id in enumerate(inputs['filtered_data'].SampleIds)
                 if id in interest_ids]
    potential_otus = [Otu(vals, name, cfg['min_abundance'], i_indexes)
                      for vals, name, md
                      in inputs['filtered_data'].iterObservations()]

    pvals = [otu.pval for otu in potential_otus]
    for otu, corrected_pval in zip(potential_otus,
                                   correct_pvalues(pvals, cfg['p_val_adj'])):
        otu.corrected_pval = corrected_pval
    # Filter down to the core
    return [otu for otu in potential_otus
            if (otu.corrected_pval <= cfg['max_p'] and
                otu.i_presence_frac >= cfg['min_frac'] and
                otu.o_presence_frac <= cfg['max_out_presence'])]


def format_results(res, cfg):
    """Format the result data as a tsv
    """
    header = '\t'.join(['OTU', 'Pval', (cfg['p_val_adj'] + ' Corrected Pval'),
                        'Interest Group Presence', 'Out Group Presence'])
    # logging.info("Results for configuration: " + cfg['name'])
    return '\n'.join([header] + [str(otu) for otu in list(sorted(res))])
