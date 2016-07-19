import webapp2
import sys
from time import localtime, strftime

from email_results import send_error_as_email
from process_data import RunPipeline
import run_config


class MainPage(webapp2.RequestHandler):
    def get(self):
        page = open('index.html', 'r')
        MAIN_PAGE_HTML = page.read()
        page.close()
        upload_url = '/'
        self.response.write(MAIN_PAGE_HTML.format(upload_url))

    def post(self):
        name = self.request.get('name')
        data = self.request.get('datafile')
        mapping_file = self.request.get('groupfile').split('\n')

        factor, group = self.request.get('group').split(':')
        NTIMES = int(self.request.get('random'))
        timestamp = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        p_val_adj = self.request.get('pvaladjmethod')
        include_out = bool(self.request.get('include_out'))
        min_abundance = float(self.request.get('min_abundance'))
        random_opt = self.request.get('random_opt')

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\n' +
                      '# of randomizations: %s\n\n\n')
                     % (factor, group, p_val_adj, NTIMES))

        params = {
            'name': name,
            'factor': factor,
            'group': group,
            'p_val_adj': p_val_adj,
            'ntimes': str(NTIMES),
            'to_email': to_email,
            'timestamp': timestamp,
            'user_args': user_args,
            'include_out': include_out,
            'random_opt': random_opt,
        }

        inputs = {
            'mapping_file': mapping_file,
            'data': data,
        }

        errors_list, mapping_dict, out_group = validate_inputs(params, inputs)

        params['out_group'] = out_group
        inputs['mapping_dict'] = mapping_dict

        params['run_cfgs'] = [
            {
                'factor': factor,
                'group': group,
                'out_group': out_group,
                'min_abundance': min_abundance,
                'name': 'interest'
            }
        ]

        if include_out:
            params['run_cfgs'].append({
                'factor': factor,
                'group': out_group,
                'out_group': group,
                'min_abundance': min_abundance,
                'name': 'out'
                }
            )

        if len(errors_list) > 0:
            send_error_as_email(timestamp, user_args, '\n'.join(errors_list),
                                to_email)
            # TODO: Return form with error marked
            sys.exit(0)

        pipeline = RunPipeline(params, inputs)
        pipeline.start()

        self.redirect('/')


def validate_inputs(params, inputs):
    mapping_file = inputs['mapping_file']
    factor = params['factor']
    group = params['group']

    labels = mapping_file[0].strip().strip('#').split(run_config.DELIM)
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
    mapping_dict = get_categ_samples_dict(mapping_file, indx_categ,
                                          indx_sampleid)
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

    if params['random_opt'] not in ['row_wise', 'column_wise',
                                    'otu_label', 'samp_annot']:
        errors_list.append('Randomization option not valid')
    return (errors_list, mapping_dict, out_group)


def get_categ_samples_dict(mapping_info_list, index_categ, indx_sampleid):
    local_dict = dict()
    for l in mapping_info_list:
        if not l or l.strip()[0] == '#':
            continue
        key, val = map(l.strip().split(run_config.DELIM).__getitem__,
                       [index_categ, indx_sampleid])
        if key in local_dict:
            local_dict[key].append(val)
        else:
            local_dict[key] = [val]
    return local_dict
