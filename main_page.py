import webapp2

from time import localtime, strftime

import string
import random
import sys

from google.appengine.api import taskqueue
from google.appengine.api.taskqueue import TaskRetryOptions

from storage import init_storage
from email_results import send_results_as_email


class MainPage(webapp2.RequestHandler):
    def get(self):
        page = open('index.html', 'r')
        MAIN_PAGE_HTML = page.read()
        page.close()
        upload_url = '/'
        self.response.write(MAIN_PAGE_HTML.format(upload_url))

    def post(self):
        otu_table_biom = self.request.get('datafile')
        group_info_list = self.request.get('groupfile').split('\n')

        DELIM = '\t'
        factor, group = self.request.get('group').split(':')
        NTIMES = int(self.request.get('random'))
        datetime = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        # added a random alpha numeric to avoid conflict with another
        # (simultaneous) user request.
        OUTPFILE = self.request.get('content') + datetime + id_generator()
        key = OUTPFILE + '~~~~' + factor + '~~~~' + group

        p_val_adj = self.request.get('pvaladjmethod')

        to_email = self.request.get('email')

        errors_list, categ_samples_dict, out_group = \
            validate_inputs(factor, group, group_info_list, DELIM, int(NTIMES),
                            OUTPFILE)

        if len(errors_list) > 0:
            user_args = (('You selected the following parameters:' +
                          '\nFactor: %s\nGroup: %s\n' +
                          'Pval correction: %s\n' +
                          '# of randomizations: %s\n\n\n')
                         % (factor, group, p_val_adj, NTIMES))
            send_results_as_email(key, user_args,
                                  '\n'.join(errors_list), to_email)
            sys.exit(0)

        # make a dict and insert in ndb as a json property
        local_string_hack_dict = {
            'otu_table_biom_key': key,
            'group_info_list': group_info_list,
            'factor': factor,
            'group': group,
            'out_group': out_group,
            'p_val_adj': p_val_adj,
            'delim': DELIM,
            'ntimes': str(NTIMES),
            'outpfile': OUTPFILE,
            'to_email': to_email,
            'categ_samples_dict': categ_samples_dict,
        }

        init_storage(key, otu_table_biom, local_string_hack_dict)

        numb_tasks = int(NTIMES)/50
        # try to break the randomizations into n tasks of 50 randomizations
        # each and then run them one by one
        # finally run the significance calculation
        # http://stackoverflow.com/questions/4224564/calling-a-script-after-tasks-queue-is-empty
        taskqueue.add(url='/process_data',
                      params={'otu_table_biom_key': key,
                              'mode': 'true'},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        taskqueue.add(url='/process_data',
                      params={'otu_table_biom_key': key,
                              'mode': 'out'},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        for i in range(numb_tasks):
            taskqueue.add(url="/process_data",
                          params={'otu_table_biom_key': key,
                                  'mode': 'random'},
                          retry_options=TaskRetryOptions(task_retry_limit=0,
                                                         task_age_limit=1),
                          countdown=1)

        taskqueue.add(url="/process_results",
                      params={'otu_table_biom_key': key,
                              'numb_tasks': numb_tasks},
                      retry_options=TaskRetryOptions(task_retry_limit=0,
                                                     task_age_limit=1),
                      countdown=1)
        self.redirect('/')


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def validate_inputs(factor, group, mapping_info_list, DELIM, NTIMES,
                    OUTPFILE):
    # find index of SampleID and category to be summarized. e.g. swg or non-swg
    labels = mapping_info_list[0].split(DELIM)
    indx_sampleid = indx_categ = ''

    errors_list = list()

    if '#SampleID' in labels:
        indx_sampleid = labels.index('#SampleID')
    else:
        errors_list.append("' not in the headers of the sample <-> " +
                           "group info file")
    if factor in labels:
        indx_categ = labels.index(factor)
    else:
        errors_list.append("'%s' not in the headers of the sample <-> " +
                           "group info file" % factor)
    if int(NTIMES) % 50 != 0:
        errors_list.append("Number of randomizations requested is not " +
                           "multiple of 50. Kindly rerun job")
    categ_samples_dict = get_categ_samples_dict(mapping_info_list, DELIM,
                                                indx_categ, indx_sampleid)
    groups = categ_samples_dict.keys()
    if (len(groups) != 2):
        errors_list.append('\nERROR: Following code divides samples ' +
                           'in >TWO groups. Change the mapping file ' +
                           'to only have two groups (e.g. A vs D)\n')
    if group not in groups:
        errors_list.append('Interest group is not present in groupfile')
        out_group = ''
    else:
        groups.remove(group)
        out_group = groups[0]
    return (errors_list, categ_samples_dict, out_group)


def get_categ_samples_dict(mapping_info_list, DELIM, index_categ,
                           indx_sampleid):
    local_dict = dict()
    for l in mapping_info_list:
        if not l or l.strip()[0] == '#':
            continue
        key, val = map(l.strip().split(DELIM).__getitem__,
                       [index_categ, indx_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict
