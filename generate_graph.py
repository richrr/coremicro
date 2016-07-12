from biom.parse import parse_biom_table
import StringIO
import matplotlib.pyplot as plt
import numpy


def generate_graph(params, inputs, results):
    mapping = inputs['mapping_dict']
    data = parse_biom_table(inputs['data'])
    attachments = list()
    cfg_to_group = {cfg['name']: cfg['group'] for cfg in params['run_cfgs']}
    for cfg in results:
        for frac in results[cfg]:
            otus = [res['otu'] for res in results[cfg][frac]]
            group = cfg_to_group[cfg]
            # Filter down to core otus
            core = data.filterObservations(lambda values, id, md:
                                           md['taxonomy'] in otus)
            interest_core = core.filterSamples(lambda values, id, md:
                                               id in mapping[group])
            interest_core_samples = len(interest_core.SampleIds)
            interest_averages = [float(total) / interest_core_samples
                                 for total in interest_core.sum(
                                         axis='observation')]
            interest_frequencies = interest_core.transformSamples(
                lambda l, id, md: numpy.array([v > 0 for v in l])
            ).sum(axis='observation')
            out_core = core.filterSamples(lambda values, id, md:
                                          id not in mapping[group])
            out_core_samples = len(out_core.SampleIds)
            out_averages = [float(total) / out_core_samples
                            for total in out_core.sum(axis='observation')]
            out_frequencies = out_core.transformSamples(
                lambda l, id, md: numpy.array([v > 0 for v in l])
            ).sum(axis='observation')

            # OTUs in the same order as the averages
            ordered_otus = [observation[2]['taxonomy']
                            for observation in core.iterObservations()]
            ordered_ids = [observation[1]
                           for observation in core.iterObservations()]

            # Sort everything by decreasing frequency in interest group
            interest_averages, interest_frequencies, out_averages, \
                out_frequencies, ordered_ids, ordered_otus \
                = zip(*reversed(sorted(zip(
                    interest_averages, interest_frequencies,
                    out_averages, out_frequencies,
                    ordered_ids, ordered_otus))))

            width = 0.35
            ind = [i + width/2 for i in range(len(otus))]

            interest = plt.bar(ind, interest_averages, width, color='r')
            out = plt.bar([i + width for i in ind], out_averages, width,
                          color='y')
            plt.ylabel('Average Abundance')
            plt.xlabel('Sample ID')
            plt.xticks([i + width for i in ind], ordered_ids)
            config = [c for c in params['run_cfgs'] if c['name'] == cfg][0]
            plt.legend((interest[0], out[0]), (config['group'],
                                               config['out_group']))
            plt.title('Abundance of Core Microbes at %s%% Threshold' % frac)

            out = StringIO.StringIO()
            plt.savefig(out, format='svg')
            attachments.append(('%s_plot_%s_%s.svg' % (cfg, frac,
                                                       params['name']),
                                out.getvalue()))

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
            plt.clf()
    return attachments
