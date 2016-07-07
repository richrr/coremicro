import numpy.random
import itertools


def randomize(mapping, data, random_opt):
    # TODO: validate that random_opt is a valid option
    if random_opt == 'row_wise':
        return row_wise(mapping, data)
    elif random_opt == 'column_wise':
        return column_wise(mapping, data)
    elif random_opt == 'otu_label':
        return otu_label(mapping, data)
    elif random_opt == 'samp_annot':
        return samp_annot(mapping, data)


def row_wise(mapping, data):
    row_randomized_data = data.copy()
    row_randomized_data.transformObservations(
        lambda data, obs_id, metadata: shuffle_list(data))
    return mapping, row_randomized_data


def column_wise(mapping, data):
    column_randomized_data = data.copy()
    column_randomized_data.transformSamples(
        lambda data, samp_id, metadata: shuffle_list(data))
    return mapping, column_randomized_data


def otu_label(mapping, data):
    ids = []
    metadata = []
    for observation in data.iterObservations():
        obs_value, obs_id, obs_metadata = observation
        ids.append(obs_id)
        metadata.append(obs_metadata)
    relabled_data = data.copy()
    relabled_data.addObservationMetadata(dict(zip(shuffle_list(ids),
                                                  metadata)))
    return mapping, relabled_data


def samp_annot(mapping, data):
    keys = mapping.keys()
    values = shuffle_list(list(itertools.chain(*mapping.values())))
    first_group_len = len(mapping.values()[0])
    randomized_mapping = dict(zip(keys, [values[:first_group_len],
                                         values[first_group_len:]]))
    return randomized_mapping, data


def shuffle_list(l):
    numpy.random.shuffle(l)
    return l
