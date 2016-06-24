import collections
import logging

from mapreduce.base_handler import PipelineBase

from storage import Result_RandomDict
from make_tree import make_tree


class ProcessResultsPipeline(PipelineBase):
    def run(self, params, core, out):
        logging.info('Processing results')
        p_val_adj = params['p_val_adj']
        DELIM = params['delim']
        NTIMES = params['ntimes']
        rand_res = Result_RandomDict.get_entries(self.root_pipeline_id)
        core_rand, out_rand = combine_results(rand_res)
        logging.info('Processing core data')
        core_res = perform_sign_calc(core_rand, p_val_adj, DELIM, core,
                                     NTIMES)
        logging.info('Processing out-group data')
        out_res = perform_sign_calc(out_rand, p_val_adj, DELIM, out,
                                    NTIMES)
        tree = make_tree(core_res, out_res)
        results_string = format_results(core_res, p_val_adj)
        return [results_string, tree]


def combine_results(results):
    core = dict()
    out = dict()
    for res in results:
        res_core = res.core
        for threshold, otus in res_core.items():
            if threshold in core:
                core[threshold].update(otus)
            else:
                core[threshold] = collections.Counter(otus)
        res_out = res.out
        for threshold, otus in res_out.items():
            if threshold in out:
                out[threshold].update(otus)
            else:
                out[threshold] = collections.Counter(otus)
    return core, out


def format_results(results, p_val_adj):
    sign_results = (('Significant results:\nOTU\tFreq. in randomized ' +
                     'data\tpval=freq/times randomized\t%s corrected pval\n')
                    % p_val_adj)
    for frac in sorted(results.keys()):
        sign_results += '\n#Frac thresh %s\n' % str(frac)
        for otu in results[frac]:
            sign_results += '%s\t%s\t%s\t%s\n' % (otu['otu'], otu['freq'],
                                                  otu['pval'],
                                                  otu['corrected_pval'])
    return sign_results


def perform_sign_calc(glob_qry_entries_in_result_rand_dict,
                      p_val_adj, DELIM, true_result_frac_thresh_otus_dict,
                      NTIMES):
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
            otus_pval = freq/float(NTIMES)
            if otus_pval < p_val:
                signif_core_microb_otu_dict[o] = otus_pval

        # check if there is at least one significant entry so far:
        if len(signif_core_microb_otu_dict) == 0:
            continue            # go to next iteration
        else:
            logging.info('There are %s signif core microb without corrections'
                         % str(len(signif_core_microb_otu_dict)))

        # adjust the pvalues of significant otus for multiple testing
        new_p_vals = list()
        if p_val_adj == 'bf':
            new_p_vals = correct_pvalues_for_multiple_testing(
                signif_core_microb_otu_dict.values(), "Bonferroni")
        elif p_val_adj == 'bh':
            new_p_vals = correct_pvalues_for_multiple_testing(
                signif_core_microb_otu_dict.values(), "Benjamini-Hochberg")
        elif p_val_adj == 'none':
            new_p_vals = signif_core_microb_otu_dict.values()

        counter = 0
        for o in signif_core_microb_otu_dict.keys():
            freq = int(randomized_otus[o])
            # p value before correction (from randomized runs)
            otus_pval = freq/float(NTIMES)
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
