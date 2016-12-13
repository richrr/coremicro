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


def row_randomize_probability(data):
    """The probability that if n_interest valsues are chosen from vals the
    number present in the interest group will be as many or greater than what
    was originally found. This is calculated with a one-tailed Fisher's Exact
    Test
    """
    return (sum([nCr(data.present, need) *
                 nCr(data.total - data.present,
                     data.n_interest - need)
                 for need in xrange(data.i_present,
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
