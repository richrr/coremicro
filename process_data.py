import logging

from email_results import send_error_as_email, send_results_as_email
from read_table import read_table
import run_config
from probability import row_randomize_probability
from generate_graph import generate_graph

import pipeline


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
        pvals = [row_randomize_probability(otu_to_vals[otu],
                                           n_interest, frac,
                                           cfg['min_abundance'])
                 for otu in true[frac]]
        pvals_corrected = correct_pvalues_for_multiple_testing(
            pvals, params['p_val_adj'])
        results[frac] = [{'otu': otu,
                          'pval': pvals[i],
                          'corrected_pval': pvals_corrected[i],
                          'threshold': frac}
                         for i, otu in enumerate(true[frac])
                         if pvals_corrected[i] <= MAX_PVAL]
    return results


def core_otus(params, inputs, cfg):
    interest_presence_fracs = {
        otu: float(sum([v > cfg['min_abundance'] for v in vals])) / len(vals)
        for vals, otu, md in inputs['filtered_data'].filterSamples(
                lambda values, id, md:
                id in inputs['mapping_dict'][cfg['group']]
        ).iterObservations()
    }
    total_presence_fracs = {
        otu: float(sum([v > cfg['min_abundance'] for v in vals])) / len(vals)
        for vals, otu, md in inputs['filtered_data'].iterObservations()
    }
    return {frac: [otu for otu in interest_presence_fracs.keys()
                   if (interest_presence_fracs[otu] >= frac and
                       total_presence_fracs[otu] < frac)]
            for frac in run_config.FRACS}


def format_results(res, params, cfg):
    attachments = list()
    sign_results = (('OTU\tpval\t%s ' +
                     'corrected pval\tthreshold\n')
                    % params['p_val_adj'])
    for frac in reversed(sorted(res.keys())):
        for otu in res[frac]:
            sign_results += '%s\t%s\t%s\t%s\n' % (
                otu['otu'], otu['pval'],
                otu['corrected_pval'], int(frac * 100))
    attachments.append(('%s_results_%s.tsv' % (cfg['name'], params['name']),
                        sign_results))
    return attachments


def correct_pvalues_for_multiple_testing(pvalues, correction_type):
    n = len(pvalues)
    if correction_type == 'bf':  # Bonferroni
        new_pvalues = [n * p for p in pvalues]
    elif correction_type == 'bf-h':  # Bonferroni-Holm
        values = sorted([(p, i) for i, p in enumerate(pvalues)])
        new_pvalues = [None] * n
        for rank, (p, i) in enumerate(values):
            new_pvalues[i] = (n - rank) * p
    elif correction_type == 'b-h':  # Benjamini-Hochberg
        values = list(reversed(sorted(
            [(p, i) for i, p in enumerate(pvalues)]
        )))
        adjusted_vals = [float(n)/(n - i) * p
                         for i, (p, index) in enumerate(values)]
        for i in xrange(n - 1):
            if adjusted_vals[i] < adjusted_vals[i + 1]:
                adjusted_vals[i + 1] = adjusted_vals[i]
        new_pvalues = [None] * n
        for i, (p, index) in enumerate(values):
            new_pvalues[index] = adjusted_vals[i]
    elif correction_type == 'none':
        new_pvalues = pvalues[:]
    else:
        logging.warn('Invalid correction type of "%s" given to correct_pvalues'
                     % correction_type)
    return new_pvalues
