def correct_pvalues(pvalues, correction_type):
    """Applies the given correction method to the given pvalues to correct for
    multiple testing errors"""
    if correction_type == 'bf':
        return bonferroni_correct(pvalues)
    elif correction_type == 'bf-h':
        return bonferroni_holm_correct(pvalues)
    elif correction_type == 'b-h':
        return benjamini_hotchberg_correct(pvalues)
    else:
        return pvalues[:]


def bonferroni_correct(pvalues):
    return [len(pvalues) * p for p in pvalues]


def bonferroni_holm_correct(pvalues):
    n = len(pvalues)
    values = sorted([(p, i) for i, p in enumerate(pvalues)])
    new_pvalues = [None] * n
    for rank, (p, i) in enumerate(values):
        new_pvalues[i] = (n - rank) * p
    return new_pvalues


def benjamini_hotchberg_correct(pvalues):
    n = len(pvalues)
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
    return new_pvalues
