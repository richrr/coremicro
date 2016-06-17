import webapp2

from time import localtime, strftime

import sys

from email_results import send_error_as_email
from process_data import RunPipeline


class MainPage(webapp2.RequestHandler):
    def get(self):
        page = open('index.html', 'r')
        MAIN_PAGE_HTML = page.read()
        page.close()
        upload_url = '/'
        self.response.write(MAIN_PAGE_HTML.format(upload_url))

    def post(self):
        data = self.request.get('datafile')
        mapping_file = self.request.get('groupfile').split('\n')

        DELIM = '\t'
        factor, group = self.request.get('group').split(':')
        NTIMES = int(self.request.get('random'))
        timestamp = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        p_val_adj = self.request.get('pvaladjmethod')

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\n' +
                      '# of randomizations: %s\n\n\n')
                     % (factor, group, p_val_adj, NTIMES))

        params = {
            'mapping_file': mapping_file,
            'factor': factor,
            'group': group,
            'p_val_adj': p_val_adj,
            'delim': DELIM,
            'ntimes': str(NTIMES),
            'to_email': to_email,
            'data': data,
            'timestamp': timestamp,
            'user_args': user_args,
        }

        errors_list, mapping_dict, out_group = validate_inputs(params)

        params['out_group'] = out_group
        params['mapping_dict'] = mapping_dict

        if len(errors_list) > 0:
            send_error_as_email(timestamp, user_args, '\n'.join(errors_list),
                                to_email)
            # TODO: Return form with error marked
            sys.exit(0)

        pipeline = RunPipeline(params)
        pipeline.start()

        self.redirect('/')


def validate_inputs(params):
    mapping_file = params['mapping_file']
    DELIM = params['delim']
    factor = params['factor']
    group = params['group']

    labels = mapping_file[0].strip().strip('#').split(DELIM)
    indx_sampleid = indx_categ = 0

    errors_list = list()

    if 'SampleID' in labels:
        indx_sampleid = labels.index('SampleID')
    else:
        errors_list.append('"SampleID" not in the headers of the sample ' +
                           '<-> group info file')
    if factor in labels:
        indx_categ = labels.index(factor)
    else:
        errors_list.append(('"%s" not in the headers of the sample <-> ' +
                            'group info file') % factor)
    mapping_dict = get_categ_samples_dict(mapping_file, DELIM,
                                          indx_categ, indx_sampleid)
    groups = mapping_dict.keys()
    if (len(groups) != 2):
        errors_list.append('Following code divides samples ' +
                           'in >TWO groups. Change the mapping file ' +
                           'to only have two groups (e.g. A vs D)\n')
    if group not in groups:
        errors_list.append('Interest group is not present in groupfile')
        out_group = ''
    else:
        groups.remove(group)
        out_group = groups[0]
    return (errors_list, mapping_dict, out_group)


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
