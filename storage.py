from google.appengine.ext import ndb
import logging


# for every random dict entry, it has different thresholds
# the actual key is automatically generated
class Results(ndb.Model):
    res = ndb.JsonProperty()

    @classmethod
    def add_entry(cls, root_id, run_id, res):
        logging.info('Writing results from pipeline %s', run_id)
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
