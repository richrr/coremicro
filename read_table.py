from biom.parse import parse_biom_table
from biom.table import table_factory, SparseOTUTable


def read_table(table_file):
    parsed_table = parse_biom_table(table_file)
    sample_ids = parsed_table.SampleIds
    samples = len(sample_ids)
    otu_data = dict()

    for vals, id, md in parsed_table.iterObservations():
        otu = md['taxonomy']
        if otu_data.get(otu):
            otu_data[otu] = [vals[i] + otu_data[otu][i]
                             for i in xrange(samples)]
        else:
            otu_data[otu] = vals

    observation_ids = otu_data.keys()
    data = [otu_data[o] for o in observation_ids]
    return table_factory(data, sample_ids, observation_ids,
                         constructor=SparseOTUTable)
