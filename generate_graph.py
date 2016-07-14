from biom.parse import parse_biom_table
import StringIO
import os
import numpy
import logging

# matplotlib can't be run on the development server
IS_NOT_DEVELOPMENT = os.getenv('SERVER_SOFTWARE',
                               '').startswith('Google App Engine/')
if IS_NOT_DEVELOPMENT:
    import matplotlib.pyplot as plt


def generate_graph(params, inputs, results):
    mapping = inputs['mapping_dict']
    data = parse_biom_table(inputs['data'])
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
                                           md['taxonomy'] in otus)
            interest_core = core.filterSamples(lambda values, id, md:
                                               id in mapping[group])
            out_core = core.filterSamples(lambda values, id, md:
                                          id not in mapping[group])
            interest_core_samples = len(interest_core.SampleIds)
            out_core_samples = len(out_core.SampleIds)

            ordered_otus = []
            ordered_ids = []

            interest_averages = []
            interest_frequencies = []
            interest_errors = []
            for vals, id, md in interest_core.iterObservations():
                ordered_otus.append(md['taxonomy'])
                ordered_ids.append(id)
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
                ordered_ids, ordered_otus \
                = zip(*reversed(sorted(zip(
                    interest_averages, interest_frequencies, interest_errors,
                    out_averages, out_frequencies, out_errors,
                    ordered_ids, ordered_otus))))

            width = 0.35
            ind = [i + width/2 for i in range(len(otus))]
            if IS_NOT_DEVELOPMENT:
                interest = plt.bar(ind, interest_averages, width, color='r',
                                   yerr=interest_errors)
                out = plt.bar([i + width for i in ind], out_averages, width,
                              color='y', yerr=out_errors)
                plt.ylabel('Average Abundance')
                plt.xlabel('Sample ID')
                plt.xticks([i + width for i in ind], ordered_ids)
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
            else:
                logging.warn('Not generating graphs because matplotlib ' +
                             'does not work in development')

            ref_text = 'ID\tInterest Frequency\tOut Frequency\tOTU\n'
            for i in range(len(otus)):
                ref_text += '%s\t%d of %d\t%d of %d\t%s\n' % (
                    ordered_ids[i],
                    interest_frequencies[i],
                    interest_core_samples,
                    out_frequencies[i],
                    out_core_samples,
                    ordered_otus[i])

            attachments.append(('%s_plot_labels_%s_%s.tsv' %
                                (cfg, frac, params['name']),
                                ref_text))
    return attachments


def standard_error(a):
    return numpy.std(a, ddof=1)/numpy.sqrt(len(a))
