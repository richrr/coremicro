from biom.parse import parse_biom_table


def read_table(table_file):
    return parse_biom_table(
        table_file
    ).collapseObservationsByMetadata(
        lambda md: md['taxonomy'],
        norm=False,
        min_group_size=1)
