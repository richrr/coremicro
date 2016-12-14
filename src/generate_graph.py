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


import StringIO
import numpy
import logging
from collections import namedtuple
import base64

from parse_inputs import samples
import run_config

# matplotlib can't be run on the development server
if run_config.IS_PRODUCTION:
    import matplotlib.pyplot as plt

# Namedtuple to hold statistics calculated for otus
Stats = namedtuple('Stats', ['otu', 'i_average', 'i_frequency', 'i_error',
                             'o_average', 'o_frequency', 'o_error'])


def generate_graph(inputs, cfg, results):
    attachments = list()
    otus = [res['otu'] for res in results]
    if len(otus) == 0:
        return attachments

    stats, i_samples, o_samples = get_stats(inputs, otus,
                                            cfg['group'],
                                            cfg['out_group'],
                                            cfg['min_abundance'])

    # Sort everything by decreasing average presence in interest group
    stats = list(reversed(sorted(
        stats, cmp=lambda a, b: cmp(a.i_average, b.i_average)
    )))

    if run_config.IS_PRODUCTION:
        attachments.append({
            'Content-Type': 'image/svg+xml',
            'Filename': '%s_plot_%s.svg' % (cfg['name'], cfg['run_name']),
            'content': base64.b64encode(make_graph(stats,
                                                   cfg['group_name'],
                                                   cfg['out_group_name']))
        })

    ref_text = 'ID\tInterest Frequency\tOut Frequency\tOTU\n'
    for i in range(len(otus)):
        ref_text += '%s\t%d of %d\t%d of %d\t%s\n' % (
            i,
            stats[i].i_frequency,
            i_samples,
            stats[i].o_frequency,
            o_samples,
            stats[i].otu)

    attachments.append({
        'Content-Type': 'text/tab-separated-values',
        'Filename': '%s_plot_labels_%s.tsv' % (cfg['name'], cfg['run_name']),
        'content': base64.b64encode(ref_text)
    })
    if not run_config.IS_PRODUCTION:
        logging.warn('Graphs not generated because in development mode')
    return attachments


def get_stats(inputs, otus, i_group, o_group, min_abundance):
    core = inputs['filtered_data'].filterObservations(
        lambda values, id, md: id in otus
    )
    interest = core.filterSamples(
        lambda values, id, md: id in samples(inputs['mapping_dict'], i_group)
    )
    out = core.filterSamples(
        lambda values, id, md: id in samples(inputs['mapping_dict'], o_group)
    )

    res = list()
    for ((i_vals, i_otu, i_md),
         (o_vals, o_id, o_md)) in zip(interest.iterObservations(),
                                      out.iterObservations()):
        res.append(Stats(i_otu,
                         numpy.mean(i_vals),
                         sum([v > min_abundance for v in i_vals]),
                         standard_error(i_vals),
                         numpy.mean(o_vals),
                         sum([v > min_abundance for v in o_vals]),
                         standard_error(o_vals)))

    i_samples = len(interest.SampleIds)
    o_samples = len(out.SampleIds)
    return res, i_samples, o_samples


def standard_error(a):
    return numpy.std(a, ddof=1)/numpy.sqrt(len(a))


def make_graph(stats, i_group_name, o_group_name):
    width = 0.35
    ind = [i + width/2 for i in range(len(stats))]
    interest = plt.bar(ind, [s.i_average for s in stats],
                       width, color='r',
                       yerr=[s.i_error for s in stats])
    out = plt.bar([i + width for i in ind],
                  [s.o_average for o in stats], width,
                  color='y', yerr=[s.o_error for s in stats])
    plt.ylabel('Average Abundance')
    plt.xlabel('Sample ID')
    plt.xticks([i + width for i in ind], range(len(stats)))
    plt.legend((interest[0], out[0]), (i_group_name, o_group_name))
    plt.title('Abundance of Core Microbes')
    out = StringIO.StringIO()
    plt.savefig(out, format='svg')
    plt.clf()
    return out.getvalue()
