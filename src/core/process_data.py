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

from correct_p_values import correct_pvalues
from probability import row_randomize_probability


def process(inputs, cfg):
    """Finds the core OTUs"""
    i_indexes = [i for i, id in enumerate(inputs['filtered_data'].SampleIds)
                 if id in [otu
                           for g in cfg['group']
                           for otu in inputs['mapping_dict'][g]]]
    potential_otus = list()
    for vals, otu, md in inputs['filtered_data'].iterObservations():
        otu_summary = make_otu_sumary(
            otu,
            total=len(vals),
            present=len(filter_min_val(vals, cfg['min_abundance'])),
            interest=len(i_indexes),
            i_present=len(filter_min_val(filter_indexes(vals, i_indexes),
                                         cfg['min_abundance'])))
        otu_summary['pval'] = row_randomize_probability(otu_summary)
        potential_otus.append(otu_summary)

    pvals = [otu['pval'] for otu in potential_otus]
    for otu, corrected_pval in zip(potential_otus,
                                   correct_pvalues(pvals, cfg['p_val_adj'])):
        otu['corrected_pval'] = corrected_pval
    # Filter down to the core
    return [otu for otu in potential_otus
            if (otu['corrected_pval'] <= cfg['max_p'] and
                otu['presence'] >= cfg['min_frac'] and
                otu['out_presence'] <= cfg['max_out_presence'])]


def make_otu_sumary(otu, total, present, interest, i_present):
    return {
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
        'presence': i_present / float(interest),
        'out_presence': (present - i_present) / float(total - interest)
    }


def filter_indexes(vals, indexes):
    """Returns a list of the values in vals that are at one of the
    indices in indexes"""
    return [v for i, v in enumerate(vals) if i in indexes]


def filter_min_val(vals, min_val):
    """Returns a list of the values in vals that are at or above the
    specified minimum abundance"""
    return [v for v in vals if v > min_val]


def format_results(res, cfg):
    """Format the result data as a tsv
    """
    def cmp_otu_results(a, b):
        """Compares results for sorting. Sorts first by descending presence,
        then by ascending corrected p-value, and finally by OTU name"""
        return (cmp(b['presence'], a['presence']) or
                cmp(a['corrected_pval'], b['corrected_pval']) or
                cmp(a['otu'], b['otu']))

    sign_results = (('OTU\tpval\t%s Corrected Pval\t' +
                     'Presence\tOut Group Presence\n')
                    % cfg['p_val_adj'])
    for otu in list(sorted(
            res, cmp=cmp_otu_results)):
        sign_results += '%s\t%s\t%s\t%s\t%s\n' % (
            otu['otu'], otu['pval'],
            otu['corrected_pval'], otu['presence'], otu['out_presence'])
    # logging.info("Results for configuration: " + cfg['name'])
    logging.info(sign_results)
    return sign_results
