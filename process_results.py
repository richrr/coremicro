import collections
import logging

import pipeline

from storage import Results
from generate_graph import generate_graph


class ProcessResultsPipeline(pipeline.Pipeline):
    def run(self, params, inputs, true_res):
        logging.info('Processing results')
        rand = combine_results(
            Results.get_entries(self.root_pipeline_id))
        res = {k: perform_sign_calc(v, rand[k], params)
               for k, v in true_res.iteritems()}
        attachments = list()
        attachments += format_results(res, params)
        attachments += generate_graph(params, inputs, res)
        return attachments


def get_attachments(res, graphs, params):
    attachments = [('%s_results_%s.tsv' % (k, params['name']),
                    format_results(v, params['p_val_adj']))
                   for k, v in res.iteritems()]
    attachments += [('%s_plot_%s.svg' % (k, params['name']), v)
                    for k, v in graphs.iteritems()]
    return attachments


def combine_results(results):
    combined = {k: dict() for k, v in results.get().res.iteritems()}
    for result in results:
        for cfg in result.res:
            for threshold, otus in result.res[cfg].items():
                if threshold in combined[cfg]:
                    combined[cfg][threshold].update(otus)
                else:
                    combined[cfg][threshold] = collections.Counter(otus)
    return combined


def format_results(res, params):
    attachments = list()
    for cfg in res:
        sign_results = (('OTU\tFreq. in randomized ' +
                         'data\tpval=freq/times randomized\t%s ' +
                         'corrected pval\tthreshold\n')
                        % params['p_val_adj'])
        for frac in reversed(sorted(map(int, res[cfg].keys()))):
            for otu in res[cfg][str(frac)]:
                sign_results += '%s\t%s\t%s\t%s\t%s\n' % (
                    otu['otu'], otu['freq'], otu['pval'],
                    otu['corrected_pval'], frac)
        print sign_results
        attachments.append(('%s_results_%s.tsv' % (cfg, params['name']),
                            sign_results))
    return attachments


def perform_sign_calc(true_result_frac_thresh_otus_dict,
                      glob_qry_entries_in_result_rand_dict,
                      params):
    p_val = 0.05
    results = {}
    logging.info('Performing significance calculations')
    logging.info('%d thresholds present in randomized data',
                 len(glob_qry_entries_in_result_rand_dict))
    for frac_s in glob_qry_entries_in_result_rand_dict.keys():
        results[frac_s] = []
        randomized_otus = glob_qry_entries_in_result_rand_dict[frac_s]
        logging.info('Number of OTUs from randomized dict for %s%% ' +
                     'threshold = %d',
                     frac_s, len(randomized_otus))

        # check significance
        signif_core_microb_otu_dict = collections.OrderedDict()
        for o in true_result_frac_thresh_otus_dict[str(frac_s)]:
            freq = 0.1
            # the else for this means it was not observed even once in the
            # randomized data!
            if o in randomized_otus:
                freq = int(randomized_otus[o])
            otus_pval = freq/float(params['ntimes'])
            if otus_pval < p_val:
                signif_core_microb_otu_dict[o] = otus_pval

        # check if there is at least one significant entry so far:
        if len(signif_core_microb_otu_dict) == 0:
            logging.info('There are no signif core microb without corrections')
            continue            # go to next iteration
        else:
            logging.info('There are %s signif core microb without corrections'
                         % str(len(signif_core_microb_otu_dict)))

        # adjust the pvalues of significant otus for multiple testing
        new_p_vals = list()
        if params['p_val_adj'] == 'bf':
            new_p_vals = correct_pvalues_for_multiple_testing(
                signif_core_microb_otu_dict.values(), "Bonferroni")
        elif params['p_val_adj'] == 'bh':
            new_p_vals = correct_pvalues_for_multiple_testing(
                signif_core_microb_otu_dict.values(), "Benjamini-Hochberg")
        elif params['p_val_adj'] == 'none':
            new_p_vals = signif_core_microb_otu_dict.values()

        counter = 0
        for o in signif_core_microb_otu_dict.keys():
            freq = int(randomized_otus[o])
            # p value before correction (from randomized runs)
            otus_pval = freq/float(params['ntimes'])
            new_p_v = new_p_vals[counter]  # p value after correction
            if new_p_v < p_val:
                results[frac_s].append({
                    'otu': o,
                    'freq': freq,
                    'pval': otus_pval,
                    'corrected_pval': new_p_v,
                    'threshold': frac_s,
                })
            counter += 1
    return results


# http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
# the pvalues do not have to be sorted, their order is maintained in the
# results
def correct_pvalues_for_multiple_testing(pvalues,
                                         correction_type='Benjamini-Hochberg'):
    """
    consistent with R - print correct_pvalues_for_multiple_
    testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071,
             0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    elif correction_type == "None":
        new_pvalues = 1 * pvalues
    return new_pvalues
