from biom.parse import parse_biom_table
from biom.table import table_factory, SparseOTUTable

import run_config


def read_table(table_file):
    parsed_table = parse_biom_table(table_file)
    sample_ids = parsed_table.SampleIds
    otu_data = dict()

    for vals, id, md in parsed_table.iterObservations():
        otu = md['taxonomy']
        if otu_data.get(otu) is not None:
            otu_data[otu] = vals + otu_data[otu]  # both are numpy arrays
        else:
            otu_data[otu] = vals

    observation_ids = otu_data.keys()
    data = [otu_data[o] for o in observation_ids]
    return table_factory(data, sample_ids, observation_ids,
                         constructor=SparseOTUTable)


def get_categ_samples_dict(mapping_info_list, factor):
    labels = mapping_info_list[0].strip().strip('#').split(run_config.DELIM)
    index_sampleid = labels.index('SampleID')
    index_categ = labels.index(factor)
    local_dict = dict()
    for l in mapping_info_list:
        if not l or l.strip()[0] == '#':
            continue
        key, val = map(l.strip().split(run_config.DELIM).__getitem__,
                       [index_categ, index_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict


def samples(mapping_dict, group):
        return [otu for g in group for otu in mapping_dict[g]]
