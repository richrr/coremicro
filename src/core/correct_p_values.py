# Copyright 2016 Richard Rodrigues, Nyle Rodgers, Mark Williams, Virginia Tech
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
