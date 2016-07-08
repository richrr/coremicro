from biom.parse import parse_biom_table
import StringIO
import numpy
import matplotlib.pyplot as plt


def generate_graph(params, inputs, results):
    mapping = inputs['mapping_dict']
    data = parse_biom_table(inputs['data'])
    all_otus = {k: get_otus_list(results[k]) for k in results}
    attachments = list()
    cfg_to_group = {cfg['name']: cfg['group'] for cfg in params['run_cfgs']}
    for cfg in all_otus:
        otus = all_otus[cfg]
        group = cfg_to_group[cfg]
        # Filter down to core otus
        interest = data.filterObservations(lambda values, id, md:
                                           md['taxonomy'] in otus)
        core_interest = interest.filterSamples(lambda values, id, md:
                                               id in mapping[group])
        core_interest_samples = len(core_interest.SampleIds)
        core_averages = [float(total) / core_interest_samples
                         for total in core_interest.sum(axis='observation')]
        core_frequencies = core_interest.reduce(lambda s, v: s + (v > 0),
                                                axis='observation')
        out_interest = interest.filterSamples(lambda values, id, md:
                                              id not in mapping[group])
        out_interest_samples = len(out_interest.SampleIds)
        out_averages = [float(total) / out_interest_samples
                        for total in out_interest.sum(axis='observation')]
        out_frequencies = out_interest.reduce(lambda s, v: s + (v > 0),
                                              axis='observation')

        # OTUs in the same order as the averages
        ordered_otus = [observation[2]['taxonomy']
                        for observation in interest.iterObservations()]
        ordered_ids = [observation[1]
                       for observation in interest.iterObservations()]

        width = 0.35
        ind = range(len(otus))

        fig, ax = plt.subplots()

        core = ax.bar(ind, core_averages, width, color='r')
        out = ax.bar([i + width for i in ind], out_averages, width, color='y')
        ax.set_ylabel('Average Presence')
        ax.set_xlabel('Sample ID')
        ax.set_xticks([i + width for i in ind])
        ax.set_xticklabels(ordered_ids, rotation='vertical')
        ax.legend((core[0], out[0]), ('Core', 'Out'))

        out = StringIO.StringIO()
        plt.savefig(out, format='svg')
        attachments.append(('%s_plot_%s.svg' % (cfg, params['name']),
                            out.getvalue()))

        ref_text = 'ID\tCore Frequency\tOut Frequency\tOTU\n'
        for i in range(len(otus)):
            ref_text += '%s\t%d of %d\t%d of %d\t%s\n' % (
                ordered_ids[i],
                core_frequencies[i],
                core_interest_samples,
                out_frequencies[i],
                out_interest_samples,
                ordered_otus[i])

        attachments.append(('%s_plot_labels_%s.tsv' % (cfg, params['name']),
                            ref_text))
        plt.clf()
    return attachments


def get_otus_list(res):
    otus_list = []
    for frac in res:
        otus_list += map(lambda d: d['otu'], res[frac])
    return list(set(otus_list))
