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
import logging

from parse_inputs import summarize_otu_data
from correct_p_values import correct_pvalues
from probability import row_randomize_probability


def process(inputs, cfg):
    """Finds the core OTUs"""
    # Get summarized data
    potential_otus = summarize_otu_data(inputs['filtered_data'],
                                        inputs['mapping_dict'], cfg['group'],
                                        cfg['min_abundance'])
    # Add pval and presence to OTUs
    for otu in potential_otus:
        otu['pval'] = row_randomize_probability(otu)
        otu['presence'] = otu['i_present'] / float(otu['interest'])
    # Add corrected pval to OTUs
    corrected_pvalues = correct_pvalues(
        [otu['pval'] for otu in potential_otus], cfg['p_val_adj']
    )
    for otu, corrected_pval in zip(potential_otus, corrected_pvalues):
        otu['corrected_pval'] = corrected_pval
    # Filter down to the core
    return [otu for otu in potential_otus
            if (otu['corrected_pval'] <= cfg['max_p'] and
                otu['presence'] >= cfg['min_frac'])]


def format_results(res, cfg):
    """Format the result data as a tsv
    """
    sign_results = (('OTU\tpval\t%s ' +
                     'corrected pval\tPresence\n')
                    % cfg['p_val_adj'])
    for otu in list(sorted(
            res, cmp=cmp_otu_results)):
        sign_results += '%s\t%s\t%s\t%s\n' % (
            otu['otu'], otu['pval'],
            otu['corrected_pval'], otu['presence'])
    # logging.info("Results for configuration: " + cfg['name'])
    logging.info(sign_results)
    return sign_results


def cmp_otu_results(a, b):
    """Compares results for sorting. Sorts first by descending presence,
    then by ascending corrected p-value, and finally by OTU name"""
    res = cmp(b['presence'], a['presence'])
    res = res if res else cmp(a['corrected_pval'], b['corrected_pval'])
    res = res if res else cmp(a['otu'], b['otu'])
    return res
