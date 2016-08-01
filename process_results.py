import logging


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
