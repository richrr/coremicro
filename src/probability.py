import math
import itertools
from operator import mul


def row_randomize_probability(vals, n_interest, frac, min_abundance=0):
    """The probability that if n_interest valsues are chosen from vals at least
    frac percent will be above min_abundance
    """
    present = sum([v > min_abundance for v in vals])
    total = len(vals)
    min_needed = int(math.ceil(frac * n_interest))
    return (sum([nCr(present, need) *
                 nCr(total - present, n_interest - need)
                 for need in xrange(min_needed, min(present,
                                                    n_interest) + 1)]) /
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


def column_randomize_probability(interest_columns, frac, min_abundance=0):
    """The probability that if one value is chosen from each column at least
    frac percent will be above min_abundance. interest_columns should be a list
    of columns
    """
    n_interest = len(interest_columns)
    # A list of the probability of a chosen value being present for each column
    probs = [float(sum([v > min_abundance for v in c]))/len(c)
             for c in interest_columns]
    # For each way of choosing enough elements such that the fraction threshold
    # can be met or exceeded calculate the probability that that combination
    # will all be present and all other columns won't be present. Then sum
    # these probabilities across all such combinations and across all possible
    # number of selectable values to meet the threshold.
    return sum([
        sum([
            product([
                probs[i] if i in comb else 1-probs[i]
                for i in range(n_interest)])
            for comb in itertools.combinations(range(n_interest), r)])
        for r in range(int(math.ceil(frac * n_interest)), n_interest + 1)])


def product(l):
    """Returns the product of the given list. Similar to sum()"""
    return reduce(mul, l)


def full_table_randomize_probability(vals, n_interest, frac, min_abundance=0):
    """The probability that if n_interest values are chosen from all values in
    the table at least frac percent will be above min_abundance. vals is
    assumed to be a two dimensional list"""
    combined_vals = [v for l in vals for v in l]
    # At this point the problems are identical
    return row_randomize_probability(combined_vals, n_interest, frac,
                                     min_abundance)
