import StringIO
import numpy
import logging

from read_table import read_table
import run_config

# matplotlib can't be run on the development server
if run_config.IS_PRODUCTION:
    import matplotlib.pyplot as plt


def generate_graph(params, inputs, results):
    mapping = inputs['mapping_dict']
    data = read_table(inputs['data'])
    attachments = list()
    cfg_to_group = {cfg['name']: cfg['group'] for cfg in params['run_cfgs']}
    for cfg in results:
        for frac in results[cfg]:
            otus = [res['otu'] for res in results[cfg][frac]]
            if len(otus) == 0:
                continue
            group = cfg_to_group[cfg]
            # Filter down to core otus
            core = data.filterObservations(lambda values, id, md:
                                           id in otus)
            for val, id, md in core.iterObservations():
                print id
            interest_core = core.filterSamples(lambda values, id, md:
                                               id in mapping[group])
            out_core = core.filterSamples(lambda values, id, md:
                                          id not in mapping[group])
            interest_core_samples = len(interest_core.SampleIds)
            out_core_samples = len(out_core.SampleIds)

            ordered_otus = []

            interest_averages = []
            interest_frequencies = []
            interest_errors = []
            for vals, id, md in interest_core.iterObservations():
                ordered_otus.append(id)
                interest_averages.append(numpy.mean(vals))
                interest_frequencies.append(numpy.count_nonzero(vals))
                interest_errors.append(standard_error(vals))

            out_averages = []
            out_frequencies = []
            out_errors = []
            for vals, id, md in out_core.iterObservations():
                out_averages.append(numpy.mean(vals))
                out_frequencies.append(numpy.count_nonzero(vals))
                out_errors.append(standard_error(vals))

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
                config = [c for c in params['run_cfgs'] if c['name'] == cfg][0]
                plt.legend((interest[0], out[0]), (config['group'],
                                                   config['out_group']))
                plt.title('Abundance of Core Microbes at %s%% Threshold' %
                          frac)

                out = StringIO.StringIO()
                plt.savefig(out, format='svg')
                attachments.append(('%s_plot_%s_%s.svg' % (cfg, frac,
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
                                (cfg, frac, params['name']),
                                ref_text))
    if not run_config.IS_PRODUCTION:
        logging.warn('Graphs not generated because in development mode')
    return attachments


def standard_error(a):
    return numpy.std(a, ddof=1)/numpy.sqrt(len(a))
