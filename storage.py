from google.appengine.ext import ndb


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

    @staticmethod
    def get_parent_key(key):
        return 'fatherresults' + key

    @classmethod
    def add_entry(cls, key, results):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, otus=results).put()

    @classmethod
    def get_entries(cls, key):
        return cls.query(cls.idx == key,
                         ancestor=ndb.Key(cls, cls.get_parent_key(key)))

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

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

    @staticmethod
    def get_parent_key(key):
        return 'fatherresultstrue' + key

    @classmethod
    def add_entry(cls, key, results):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, true_results=results).put()

    @classmethod
    def get_entry(cls, key):
        result = cls.query(cls.idx == key,
                           ancestor=ndb.Key(cls, cls.get_parent_key(key)))
        if int(result.count()) == 1:
            print "Read single entry from Truedict datastore!"
        else:
            # do something useful here
            print "Single entry not found from Truedict datastore, some error!"
        return result.get()

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

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
    def add_entry(cls, key, biom, params_str):
        return cls(parent=ndb.Key(cls, cls.get_parent_key(key)),
                   idx=key, biom=biom, params_str=params_str).put()

    @classmethod
    def get_entry(cls, key):
        result = cls.query(cls.idx == key,
                           ancestor=ndb.Key(cls, cls.get_parent_key(key)))
        if int(result.count()) == 1:
            print "Read single entry from OriginalBiom datastore!"
        else:
            # do something useful here
            print ('Single entry not found from OriginalBiom datastore, ' +
                   'some error!')
        return result.get()

    @classmethod
    def add_parent(cls, key):
        return cls(id=cls.get_parent_key(key)).put()

    @classmethod
    def delete_entries(cls, key):
        parent_key = ndb.Key(cls, cls.get_parent_key(key))
        entries = cls.query(cls.idx == key,
                            ancestor=parent_key)
        for entry in entries:
            entry.key.delete()
        parent_key.delete()
