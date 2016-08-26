import logging
import pipeline

import run_config
from email_results import send_error_as_email, send_results_as_email
from parse_inputs import read_table, get_data_summary
from probability import row_randomize_probability
from generate_graph import generate_graph
from correct_p_values import correct_pvalues


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        # Must be processed here because the otu_table class can't be given to
        # a pipeline as an argument
        inputs['filtered_data'] = read_table(inputs['data'])

        attachments = list()
        for cfg in params['run_cfgs']:
            data = get_data_summary(inputs['filtered_data'],
                                    inputs['mapping_dict'],
                                    cfg['group'],
                                    cfg['min_abundance'])
            results = get_signif_otus(params, data, cfg)
            attachments += format_results(results, params, cfg)
            attachments += generate_graph(params, inputs, cfg, results)
        send_results_as_email(params, attachments)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]

        if self.was_aborted:    # There's probably a bug in the code
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(params, error)
        self.cleanup()


def get_signif_otus(params, data, cfg):
    """Takes the set of OTUs at each fractional threshold that satisfy the
    inclusion requirement that they be present at at least the threshold in the
    interest group and at less than that threshold across all samples,
    calculates the p-value for each OTU, and gives the list of OTUs found to be
    significant after correction for multiple testing"""
    results = dict()
    core = {frac: [otu for otu in data.keys()
                   if (data[otu].i_present /
                       float(data[otu].n_interest) >= frac) and
                   data[otu].present / float(data[otu].total) < frac]
            for frac in run_config.FRACS}
    for frac in core:
        logging.info('Starting frac %f' % frac)
        pvals = [row_randomize_probability(data[otu], frac)
                 for otu in core[frac]]
        pvals_corrected = correct_pvalues(pvals, params['p_val_adj'])
        results[frac] = [{'otu': otu,
                          'pval': pvals[i],
                          'corrected_pval': pvals_corrected[i],
                          'threshold': frac}
                         for i, otu in enumerate(core[frac])
                         if pvals_corrected[i] <= params['max_p']]
    return results


def format_results(res, params, cfg):
    """Format the result data as a tsv
    """
    sign_results = (('OTU\tpval\t%s ' +
                     'corrected pval\tthreshold\n')
                    % params['p_val_adj'])
    for frac in reversed(sorted(res.keys())):
        for otu in res[frac]:
            sign_results += '%s\t%s\t%s\t%s\n' % (
                otu['otu'], otu['pval'],
                otu['corrected_pval'], int(frac * 100))
    if not run_config.IS_PRODUCTION:
        print "Results for configuration: " + cfg['name']
        print sign_results
    return [('%s_results_%s.tsv' % (cfg['name'], params['name']),
            sign_results)]
