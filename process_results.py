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
    print sign_results
    return attachments


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
    if correction_type == 'bf':  # Bonferroni
        new_pvalues = n * pvalues
    elif correction_type == 'bf-h':  # Bonferroni-Holm
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == 'b-h':  # Benjamini-Hochberg
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
    elif correction_type == "none":
        new_pvalues = 1 * pvalues
    return new_pvalues
