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
from parse_inputs import get_categ_samples_dict, read_table
from operator import mul
import itertools
import math


def product(l):
    """Returns the product of the given list. Similar to sum()"""
    return reduce(mul, l)


factor = 'Plant'
group = 'Sw'
frac = 1.0
with open('../sample-sheet-qiime.txt') as f:
    mapping = get_categ_samples_dict(f.read().split('\n'), factor)

with open('../otu_table_16s.biom') as f:
    data = read_table(f.read())

interest = data.filterSamples(lambda values, id, md: id in mapping[group])

otu = 'k__Bacteria;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Cryomorphaceae;g__;s__'

n_interest = len(mapping[group])


# print 'Row Wise'
# data_otu_to_vals = {otu: vals for vals, otu, md
#                     in data.iterObservations()}
# interest_otu_to_vals = {otu: vals for vals, otu, md
#                         in interest.iterObservations()}
# print ('%d of %d in all, %d of %d in interest' %
#        (sum([v > 0 for v in data_otu_to_vals[otu]]),
#         len(data_otu_to_vals[otu]),
#         sum([v > 0 for v in interest_otu_to_vals[otu]]),
#         len(interest_otu_to_vals[otu])))
# print row_randomize_probability(data_otu_to_vals[otu],
#                                 n_interest,
#                                 frac)

# print 'Column Wise'
# interest_columns = [vals for vals, id, md in data.iterSamples()
#                     if id in mapping[group]]
# probs = [float(sum([v > 0 for v in c]))/len(c)
#          for c in interest_columns]
# print sum([
#         sum([
#             product([
#                 probs[i] if i in comb else 1-probs[i]
#                 for i in range(n_interest)])
#             for comb in itertools.combinations(range(n_interest), r)])
#         for r in range(int(math.ceil(frac * n_interest)), n_interest + 1)])
# print product(probs)
# print column_randomize_probability(interest_columns, frac)

# print 'Full Table'
# table_vals = [vals for vals, id, md in data.iterObservations()]
# print full_table_randomize_probability(table_vals, n_interest, frac)
