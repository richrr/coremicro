import webapp2

from google.appengine.ext import ndb

from time import localtime, strftime

import string
import random

from google.appengine.api import taskqueue
from google.appengine.api.taskqueue import TaskRetryOptions

from storage import Result_RandomDict, Result_TrueDict, OriginalBiom


class MainPage(webapp2.RequestHandler):
    def get(self):
        page = open('index.html', 'r')
        MAIN_PAGE_HTML = page.read()
        page.close()
        upload_url = '/'
        self.response.write(MAIN_PAGE_HTML.format(upload_url))

    def post(self):
        otu_table_biom = self.request.get('datafile')
        group_info = self.request.get('groupfile')

        DELIM = '\t'
        factor, group = self.request.get('group').split(':')
        out_group = self.request.get('out_group')
        NTIMES = int(self.request.get('random'))
        datetime = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        # added a random alpha numeric to avoid conflict with another
        # (simultaneous) user request.
        OUTPFILE = self.request.get('content') + datetime + id_generator()

        p_val_adj = self.request.get('pvaladjmethod')

        to_email = self.request.get('email')

        # this is to query all entries in this run
        ndb_custom_key_o = OUTPFILE + '~~~~' + factor + '~~~~' + group

        # make a dict and insert in ndb as a json property
        local_string_hack_dict = {
            'otu_table_biom_key': ndb_custom_key_o,
            'g_info_not_list': group_info,
            'factor': factor,
            'group': group,
            'out_group': out_group,
            'p_val_adj': p_val_adj,
            'delim': DELIM,
            'ntimes': str(NTIMES),
            'outpfile': OUTPFILE,
            'to_email': to_email,
        }
        # temporary hack to clear out the Datastore
        # http://stackoverflow.com/questions/1062540/how-to-delete-all-datastore-in-google-app-engine
        # ndb.delete_multi(RandomDict.query().fetch(keys_only=True))
        # ndb.delete_multi(ResultFile.query().fetch(keys_only=True))
        # uncomment this to clear out datastore every time
        ndb.delete_multi(Result_RandomDict.query().fetch(keys_only=True))
        ndb.delete_multi(OriginalBiom.query().fetch(keys_only=True))
        ndb.delete_multi(Result_TrueDict.query().fetch(keys_only=True))

        origb_idx = 'origbiom' + ndb_custom_key_o
        OriginalBiom(id=origb_idx).put()  # the datastore of original biom
        OriginalBiom(parent=ndb.Key(OriginalBiom, origb_idx),
                     idx=ndb_custom_key_o, biom=otu_table_biom,
                     params_str=local_string_hack_dict).put()

        fatherres_idx = 'fatherresults' + ndb_custom_key_o
        # the datastore of results from random dicts
        Result_RandomDict(id=fatherres_idx).put()

        fatherrestrue_idx = 'fatherresultstrue' + ndb_custom_key_o
        # the datastore of results from random dicts
        Result_TrueDict(id=fatherrestrue_idx).put()

        numb_tasks = int(NTIMES)/50
        # try to break the randomizations into n tasks of 50 randomizations
        # each and then run them one by one
        # finally run the significance calculation
        # http://stackoverflow.com/questions/4224564/calling-a-script-after-tasks-queue-is-empty
        taskqueue.add(url='/process_data',
                      params={'otu_table_biom_key': ndb_custom_key_o,
                              'mode': 'true'},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        taskqueue.add(url='/process_data',
                      params={'otu_table_biom_key': ndb_custom_key_o,
                              'mode': 'out'},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        for i in range(numb_tasks):
            taskqueue.add(url="/process_data",
                          params={'otu_table_biom_key': ndb_custom_key_o,
                                  'mode': 'random'},
                          retry_options=TaskRetryOptions(task_retry_limit=0,
                                                         task_age_limit=1),
                          countdown=1)

        taskqueue.add(url="/process_results",
                      params={'otu_table_biom_key': ndb_custom_key_o,
                              'numb_tasks': numb_tasks},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        self.redirect('/')


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))
