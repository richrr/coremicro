from google.appengine.ext import ndb

import logging


def init_storage(key, biom, params_str):
    '''
    Initialize storage with parents and an OriginalBiom entry with the biom
    and parameters
    '''
    OriginalBiom.add_parent(key)
    OriginalBiom.add_entry(key, biom, params_str)
    Result_RandomDict.add_parent(key)
    Result_TrueDict.add_parent(key)


def clean_storage(key):
    '''
    Clean out all storage used with the specified key
    '''
    Result_RandomDict.delete_entries(key)
    Result_TrueDict.delete_entries(key)
    OriginalBiom.delete_entries(key)


# for every random dict entry, it has different thresholds
# the actual key is automatically generated
class Result_RandomDict(ndb.Model):
    idx = ndb.StringProperty()  # the run id
    # the core otus from the (shuffled) dictionary as a json
    otus = ndb.JsonProperty()
    is_out_group = ndb.BooleanProperty(default=False)

    @staticmethod
    def get_parent_key(key):
        return 'fatherresults' + key

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

    @classmethod
    def add_entry(cls, key, results, out_group=False):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, otus=results, is_out_group=out_group).put()

    @classmethod
    def get_entries(cls, key, out_group=False):
        return cls.query(cls.idx == key, cls.is_out_group == out_group,
                         ancestor=ndb.Key(cls, cls.get_parent_key(key)))

    @classmethod
    def delete_entries(cls, key):
        parent_key = ndb.Key(cls, cls.get_parent_key(key))
        entries = cls.query(cls.idx == key,
                            ancestor=parent_key)
        for entry in entries:
            entry.key.delete()
        parent_key.delete()


class Result_TrueDict(ndb.Model):
    idx = ndb.StringProperty()  # the run id
    # the core biom from the true dictionary as a json
    true_results = ndb.JsonProperty()
    is_out_group = ndb.BooleanProperty(default=False)

    @staticmethod
    def get_parent_key(key):
        return 'fatherresultstrue' + key

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

    @classmethod
    def add_entry(cls, key, results, out_group=False):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, true_results=results, is_out_group=out_group).put()

    @classmethod
    def get_entry(cls, key, out_group=False):
        result = cls.query(cls.idx == key, cls.is_out_group == out_group,
                           ancestor=ndb.Key(cls, cls.get_parent_key(key)))
        if int(result.count()) == 1:
            logging.info('Read single entry from Truedict datastore!')
        else:
            # do something useful here
            logging.warning(
                'Single entry not found from Truedict datastore, some error!')
        return result.get()

    @classmethod
    def delete_entries(cls, key):
        parent_key = ndb.Key(cls, cls.get_parent_key(key))
        entries = cls.query(cls.idx == key,
                            ancestor=parent_key)
        for entry in entries:
            entry.key.delete()
        parent_key.delete()


# the actual key is automatically generated
class OriginalBiom(ndb.Model):
    idx = ndb.StringProperty()  # the run id
    # the core biom from the original input as a json
    biom = ndb.JsonProperty()
    params_str = ndb.JsonProperty()

    @staticmethod
    def get_parent_key(key):
        return 'origbiom' + key

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

    @classmethod
    def add_entry(cls, key, biom, params_str):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, biom=biom, params_str=params_str).put()

    @classmethod
    def get_entry(cls, key):
        result = cls.query(cls.idx == key,
                           ancestor=ndb.Key(cls, cls.get_parent_key(key)))
        if int(result.count()) == 1:
            logging.info('Read single entry from OriginalBiom datastore!')
        else:
            # do something useful here
            logging.warning(
                'Single entry not found from OriginalBiom datastore, ' +
                'some error!')
        return result.get()

    @classmethod
    def delete_entries(cls, key):
        parent_key = ndb.Key(cls, cls.get_parent_key(key))
        entries = cls.query(cls.idx == key,
                            ancestor=parent_key)
        for entry in entries:
            entry.key.delete()
        parent_key.delete()

    @classmethod
    def get_params(cls, key):
        q_dict = cls.get_entry(key).to_dict()

        params = q_dict['params_str']  # this is a dictionary
        factor = params['factor']
        group = params['group']
        out_group = params['out_group']
        p_val_adj = params['p_val_adj']
        DELIM = params['delim']
        NTIMES = int(params['ntimes'])
        OUTPFILE = params['outpfile']
        to_email = params['to_email']

        otu_table_biom = q_dict['biom']
        group_info_list = params['group_info_list']
        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\n' +
                      '# of randomizations: %d\n\n\n')
                     % (factor, group, p_val_adj, NTIMES))
        categ_samples_dict = params['categ_samples_dict']

        return (user_args, to_email, p_val_adj, DELIM, NTIMES,
                otu_table_biom, group_info_list, factor, group, out_group,
                OUTPFILE, categ_samples_dict)
