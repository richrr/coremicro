from process_results import (format_results,
                             correct_pvalues_for_multiple_testing)

import logging

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

        attachments = list()
        for cfg in params['run_cfgs']:
            results = get_signif_otus(params, inputs, cfg,
                                      core_otus(params, inputs, cfg))
            attachments += format_results(results, params, cfg)
            attachments += generate_graph(params, inputs, cfg, results)
        send_results_as_email(params, attachments)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]

        if self.was_aborted:
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(params, error)
        self.cleanup()


def get_signif_otus(params, inputs, cfg, true_res):
    otu_to_vals = {otu: vals for vals, otu, md
                   in inputs['filtered_data'].iterObservations()}
    MAX_PVAL = 0.05
    results = dict()

    n_interest = len(inputs['mapping_dict'][cfg['group']])
    results[cfg['name']] = dict()
    true = true_res
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


def core_otus(params, inputs, cfg):
    interest_presence_fracs = {
        otu: float(sum([v > cfg['min_abundance'] for v in vals])) / len(vals)
        for vals, otu, md in inputs['filtered_data'].filterSamples(
                lambda values, id, md:
                id in inputs['mapping_dict'][cfg['group']]
        ).iterObservations()
    }
    out_presence_fracs = {
        otu: float(sum([v > cfg['min_abundance'] for v in vals])) / len(vals)
        for vals, otu, md in inputs['filtered_data'].filterSamples(
                lambda values, id, md:
                id in inputs['mapping_dict'][cfg['out_group']]
        ).iterObservations()
    }
    return {frac: [otu for otu in interest_presence_fracs.keys()
                   if (interest_presence_fracs[otu] >= frac and
                       out_presence_fracs[otu] < frac)]
            for frac in run_config.FRACS}
