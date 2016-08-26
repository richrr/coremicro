import math


def row_randomize_probability(data, frac):
    """The probability that if n_interest valsues are chosen from vals at least
    frac percent will be above min_abundance
    """
    min_needed = int(math.ceil(frac * data.n_interest))
    return (sum([nCr(data.present, need) *
                 nCr(data.total - data.present,
                     data.n_interest - need)
                 for need in xrange(min_needed,
                                    min(data.present,
                                        data.n_interest) + 1)]) /
            nCr(data.total, data.n_interest))


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
