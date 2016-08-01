import StringIO
import numpy
import logging

from read_table import read_table
import run_config

# matplotlib can't be run on the development server
if run_config.IS_PRODUCTION:
    import matplotlib.pyplot as plt


def generate_graph(params, inputs, cfg, results):
    attachments = list()
    for frac in results:
        otus = [res['otu'] for res in results[frac]]
        if len(otus) == 0:
            continue

        ordered_otus, interest_averages, interest_frequencies,\
            interest_errors, interest_core_samples = get_stats(
                inputs, otus, cfg['group'], cfg['min_abundance'])
        ordered_otus, out_averages, out_frequencies,\
            out_errors, out_core_samples = get_stats(
                inputs, otus, cfg['out_group'], cfg['min_abundance'])

        # Sort everything by decreasing frequency in interest group
        interest_averages, interest_frequencies, interest_errors, \
            out_averages, out_frequencies, out_errors, \
            ordered_otus \
            = zip(*reversed(sorted(zip(
                interest_averages, interest_frequencies, interest_errors,
                out_averages, out_frequencies, out_errors,
                ordered_otus))))

        width = 0.35
        ind = [i + width/2 for i in range(len(otus))]
        if run_config.IS_PRODUCTION:
            interest = plt.bar(ind, interest_averages, width, color='r',
                               yerr=interest_errors)
            out = plt.bar([i + width for i in ind], out_averages, width,
                          color='y', yerr=out_errors)
            plt.ylabel('Average Abundance')
            plt.xlabel('Sample ID')
            plt.xticks([i + width for i in ind], range(len(ordered_otus)))
            plt.legend((interest[0], out[0]), (cfg['group'],
                                               cfg['out_group']))
            plt.title('Abundance of Core Microbes at %s%% Threshold' %
                      int(frac * 100))

            out = StringIO.StringIO()
            plt.savefig(out, format='svg')
            attachments.append(('%s_plot_%s_%s.svg' % (cfg['name'],
                                                       int(frac * 100),
                                                       params['name']),
                                out.getvalue()))
            plt.clf()

        ref_text = 'ID\tInterest Frequency\tOut Frequency\tOTU\n'
        for i in range(len(otus)):
            ref_text += '%s\t%d of %d\t%d of %d\t%s\n' % (
                i,
                interest_frequencies[i],
                interest_core_samples,
                out_frequencies[i],
                out_core_samples,
                ordered_otus[i])

        attachments.append(('%s_plot_labels_%s_%s.tsv' %
                            (cfg['name'], int(frac * 100), params['name']),
                            ref_text))
    if not run_config.IS_PRODUCTION:
        logging.warn('Graphs not generated because in development mode')
    return attachments


def get_stats(inputs, otus, group, min_abundance):
    core = inputs['filtered_data'].filterObservations(
        lambda values, id, md: id in otus
    ).filterSamples(
        lambda values, id, md: id in inputs['mapping_dict'][group]
    )

    core_samples = len(core.SampleIds)

    ordered_otus = list()
    averages = list()
    frequencies = list()
    errors = list()

    for vals, id, md in core.iterObservations():
            ordered_otus.append(id)
            averages.append(numpy.mean(vals))
            frequencies.append(sum([v > min_abundance
                                    for v in vals]))
            errors.append(standard_error(vals))

    return ordered_otus, averages, frequencies, errors, core_samples


def standard_error(a):
    return numpy.std(a, ddof=1)/numpy.sqrt(len(a))
