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


import StringIO
import logging
import base64

import web_config

# matplotlib can't be run on the development server
if web_config.IS_PRODUCTION:
    import matplotlib.pyplot as plt


def samples(mapping_dict, group):
    return [otu for g in group for otu in mapping_dict[g]]


def generate_graph(inputs, cfg, otus):
    attachments = list()
    if not otus:
        return attachments

    # Sort everything by decreasing average presence in interest group
    otus = list(reversed(sorted(
        otus, cmp=lambda a, b: cmp(a.interest_mean, b.interest_mean)
    )))

    if web_config.IS_PRODUCTION:
        attachments.append({
            'Content-Type': 'image/svg+xml',
            'Filename': '%s_plot_%s.svg' % (cfg['name'], cfg['run_name']),
            'content': base64.b64encode(make_graph(otus,
                                                   cfg['group_name'],
                                                   cfg['out_group_name']))
        })

    ref_text = 'ID\tInterest Frequency\tOut Frequency\tOTU\n'
    for i in range(len(otus)):
        ref_text += '%s\t%d of %d\t%d of %d\t%s\n' % (
            i,
            otus[i].interest_present,
            otus[i].interest,
            otus[i].out_present,
            otus[i].out,
            otus[i].name)

    attachments.append({
        'Content-Type': 'text/tab-separated-values',
        'Filename': '%s_plot_labels_%s.tsv' % (cfg['name'], cfg['run_name']),
        'content': base64.b64encode(ref_text)
    })
    if not web_config.IS_PRODUCTION:
        logging.warn('Graphs not generated because in development mode')
    return attachments


def make_graph(otus, i_group_name, o_group_name):
    width = 0.35
    ind = [i + width/2 for i in range(len(otus))]
    interest = plt.bar(ind, [s.interest_mean for s in otus],
                       width, color='r',
                       yerr=[s.interest_error for s in otus])
    out = plt.bar([i + width for i in ind],
                  [s.out_mean for o in otus], width,
                  color='y', yerr=[s.out_error for s in otus])
    plt.ylabel('Average Abundance')
    plt.xlabel('Sample ID')
    plt.xticks([i + width for i in ind], range(len(otus)))
    plt.legend((interest[0], out[0]), (i_group_name, o_group_name))
    plt.title('Abundance of Core Microbes')
    out = StringIO.StringIO()
    plt.savefig(out, format='svg')
    plt.clf()
    return out.getvalue()
