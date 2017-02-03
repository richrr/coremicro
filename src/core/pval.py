# Copyright 2016, 2017 Richard Rodrigues, Nyle Rodgers, Mark Williams,
# Virginia Tech
#
# This file is part of Coremic.
#
# Coremic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Coremic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Coremic. If not, see <http://www.gnu.org/licenses/>.


def getpval(otu):
    """The probability that if n_interest valsues are chosen from vals the
    number present in the interest group will be as many or greater than what
    was originally found. This is calculated with a one-tailed Fisher's Exact
    Test
    """
    return sum(
        [nCr(otu.present, need) * nCr(otu.absent, otu.interest - need)
         for need in xrange(otu.interest_present,
                            min(otu.present, otu.interest) + 1)]
    ) / nCr(otu.total, otu.interest)


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
