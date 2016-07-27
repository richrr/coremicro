from itertools import compress

from process_results import (format_results,
                             correct_pvalues_for_multiple_testing)

import logging
import numpy

from email_results import send_error_as_email, send_results_as_email
from read_table import read_table
import run_config
from probability import row_randomize_probability
from generate_graph import generate_graph


import pipeline
import pipeline.common


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        inputs['filtered_data'] = read_table(inputs['data'])

        results = for_each_config(get_signif_otus,
                                  params, inputs,
                                  for_each_config(get_core_otus,
                                                  params, inputs))
        attachments = list()
        attachments += format_results(results, params)
        attachments += generate_graph(params, inputs, results)
        send_results_as_email(params, attachments)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]

        if self.was_aborted:
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(params, error)
        self.cleanup()


def for_each_config(f, params, inputs, *args):
    return {cfg['name']: f(params, inputs, cfg, *args)
            for cfg in params['run_cfgs']}


def get_signif_otus(params, inputs, cfg, true_res):
    otu_to_vals = {otu: vals for vals, otu, md
                   in inputs['filtered_data'].iterObservations()}
    MAX_PVAL = 0.05
    results = dict()

    n_interest = len(inputs['mapping_dict'][cfg['group']])
    results[cfg['name']] = dict()
    true = true_res[cfg['name']]
    for frac in true:
        results[frac] = list()
        pvals = [row_randomize_probability(otu_to_vals[otu],
                                           n_interest, frac,
                                           cfg['min_abundance'])
                 for otu in true[frac]]
        pvals_corrected = correct_pvalues_for_multiple_testing(
            pvals, params['p_val_adj'])
        for i, otu in enumerate(true[frac]):
            if pvals_corrected[i] <= MAX_PVAL:
                results[frac].append({
                    'otu': otu,
                    'pval': pvals[i],
                    'corrected_pval': pvals_corrected[i],
                    'threshold': frac,
                })
    return results


def get_random_results(params, inputs, cfg):
    results = dict()
    n_interest = len(inputs['mapping_dict'][cfg['group']])
    results = {frac: dict()
               for frac in run_config.FRACS}
    for vals, otu, md in inputs['filtered_data'].iterObservations():
        for frac, pval in zip(run_config.FRACS,
                              row_randomize_probability(
                                  vals, n_interest,
                                  cfg['min_abundance'])):
            results[frac][otu] = pval
    return results


def get_core_otus(params, inputs, cfg):
    mapping = inputs['mapping_dict']
    # a table of presence/absence data for just the interest group samples
    interest = inputs['filtered_data'].filterSamples(
        lambda values, id, md: id in mapping[cfg['group']])
    interest_samples = len(interest.SampleIds)
    presence_counts = interest.transformSamples(
        lambda l, id, md: numpy.array([v > cfg['min_abundance'] for v in l])
    ).sum(axis='observation')
    presence_fracs = [float(count) / interest_samples
                      for count in presence_counts]
    otus = [observation[1] for observation in interest.iterObservations()]
    return {frac:
            list(compress(otus, [presence >= frac
                                 for presence in presence_fracs]))
            for frac in run_config.FRACS}
