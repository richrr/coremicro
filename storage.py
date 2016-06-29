from google.appengine.ext import ndb


def clean_storage(key):
    '''
    Clean out all storage used with the specified key
    '''
    Result_RandomDict.delete_entries(key)


# for every random dict entry, it has different thresholds
# the actual key is automatically generated
class Result_RandomDict(ndb.Model):
    res = ndb.JsonProperty()

    @classmethod
    def add_entry(cls, root_id, run_id, res):
        return cls(parent=ndb.Key(cls, root_id),
                   id=run_id, res=res).put()

    @classmethod
    def get_entries(cls, root_id):
        return cls.query(ancestor=ndb.Key(cls, root_id))

    @classmethod
    def delete_entries(cls, root_id):
        entries = cls.get_entries(root_id)
        for entry in entries:
            entry.key.delete()
