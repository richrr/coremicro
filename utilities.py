from google.appengine.ext import ndb

from storage import OriginalBiom


# get unique elements from last column (otus)
def compile_results(otus, DELIM):
    taxon_list = list()         # this may be a json or list
    result_ = otus              # json.loads(otus)
    for l in result_:
        l = l.strip()
        contents = l.split(DELIM)
        if '#' in l or not l:
            continue
        taxon_list.append(contents[1])
    return list(set(taxon_list))


def get_required_params_from_orig_dict(otu_table_biom_o):
    origb_idx = 'origbiom' + otu_table_biom_o
    qry_entries_in_origbiom = OriginalBiom.query(
        OriginalBiom.idx == otu_table_biom_o,
        ancestor=ndb.Key(OriginalBiom, origb_idx))

    if int(qry_entries_in_origbiom.count()) == 1:
        print "Read single entry from OriginalBiom datastore!"
    else:
        # do something useful here
        print ('Single entry not found from OriginalBiom datastore, ' +
               'some error!')
    for q in qry_entries_in_origbiom:
        print 'Running 1st attempt'
        q_dict = q.to_dict()

        params = q_dict['params_str']  # this is a dictionary
        factor = params['factor']
        group = params['group']
        out_group = params['out_group']
        p_val_adj = params['p_val_adj']
        DELIM = params['delim']
        NTIMES = str(50)    # params['ntimes'])
        OUTPFILE = params['outpfile']
        to_email = params['to_email']

        otu_table_biom = q_dict['biom']
        g_info_list = params['g_info_not_list'].split('\n')
        user_args = (('You selected the following parameters:\n' +
                      'Factor: %s\nGroup: %s\nPval correction: %s')
                     % (factor, group, p_val_adj))

        return (user_args, to_email, p_val_adj, DELIM, NTIMES,
                otu_table_biom, g_info_list, factor, group, out_group,
                OUTPFILE)
