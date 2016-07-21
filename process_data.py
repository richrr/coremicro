from itertools import compress

from process_results import get_final_results, format_results

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

        true_res = run_data(inputs['mapping_dict'],
                            read_table(inputs['data']),
                            params['run_cfgs'])
        pval_res = get_random_results(params, inputs)
        results = get_final_results(true_res, pval_res, params)
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


def get_random_results(params, inputs):
    data = read_table(inputs['data'])
    results = dict()
    for cfg in params['run_cfgs']:
        n_interest = len(inputs['mapping_dict'][cfg['group']])
        results[cfg['name']] = {int(100 * frac): dict()
                                for frac in run_config.FRACS}
        for vals, otu, md in data.iterObservations():
            for frac, pval in zip(run_config.FRACS,
                                  row_randomize_probability(
                                      vals, n_interest,
                                      cfg['min_abundance'])):
                results[cfg['name']][int(frac * 100)][otu] = pval
    return results


def run_data(mapping, data, run_cfgs):
    res = dict()
    for cfg in run_cfgs:
        res[cfg['name']] = get_core(mapping, data, cfg['group'],
                                    min_abundance=cfg['min_abundance'])
    return res


def get_core(mapping, data, group, min_abundance=0):
    # a table of presence/absence data for just the interest group samples
    interest = data.filterSamples(lambda values, id, md: id in mapping[group])
    interest_samples = len(interest.SampleIds)
    presence_counts = interest.transformSamples(
        lambda l, id, md: numpy.array([v > min_abundance for v in l])
    ).sum(axis='observation')
    presence_fracs = [float(count) / interest_samples
                      for count in presence_counts]
    otus = [observation[1] for observation in interest.iterObservations()]
    return {int(frac * 100):
            list(compress(otus, [presence >= frac
                                 for presence in presence_fracs]))
            for frac in run_config.FRACS}
