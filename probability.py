import math


def row_randomize_probability(vals, n_interest, frac, min_abundance=0):
    """The probability that if n_interest valsues are chosen from vals at least
    frac percent will be above min_abundance for each frac in run_config.FRACS
    """
    present = sum([v > min_abundance for v in vals])
    total = len(vals)
    needed = int(math.ceil(frac * n_interest))
    if present < needed:
        return 0
    else:
        return (sum([nCr(present, need) *
                    nCr(total - present, n_interest - need)
                    for need in xrange(needed,
                                       min(present,
                                           n_interest) + 1)
                    if total - present >= n_interest - need]) /
                nCr(total, n_interest))


def nCr(n, r):
    """n Choose r
    Uses floating point arithemetic
    """
    assert 0 <= n
    assert 0 <= r and r <= n
    num = 1.0
    denom = 1.0
    for t in xrange(1, min(r, n - r) + 1):
        num *= n
        denom *= t
        n -= 1
    return num / denom
