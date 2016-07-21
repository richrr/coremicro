import math
import run_config


def row_randomize_probability(vals, n_interest, min_abundance=0):
    """The probability that if n_interest valsues are chosen from vals at least
    frac percent will be above min_abundance for each frac in run_config.FRACS
    """
    present_vals = sum([v > min_abundance for v in vals])
    total_vals = len(vals)
    results = []
    for frac in run_config.FRACS:
        needed_vals = int(math.ceil(frac * n_interest))
        if present_vals < needed_vals:
            results.append(0)
        else:
            results.append(nCr(present_vals, needed_vals) /
                           nCr(total_vals, n_interest))
    return results


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
