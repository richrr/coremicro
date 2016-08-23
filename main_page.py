import webapp2
from time import localtime, strftime

from email_results import send_error_as_email
from process_data import RunPipeline
from parse_inputs import get_categ_samples_dict
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

        factor = self.request.get('factor')
        group = map(lambda s: s.strip(), self.request.get('group').split(','))
        group_name = self.request.get('group_name')
        out_group_name = self.request.get('out_group_name')
        timestamp = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        randomization = self.request.get('randomization')
        p_val_adj = self.request.get('pvaladjmethod')
        include_out = bool(self.request.get('include_out'))
        min_abundance = float(self.request.get('min_abundance'))

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      '\nRandomization: %s\n' +
                      'Pval correction: %s\n\n\n')
                     % (factor, group, randomization, p_val_adj))

        params = {
            'name': name,
            'factor': factor,
            'group': group,
            'p_val_adj': p_val_adj,
            'to_email': to_email,
            'timestamp': timestamp,
            'user_args': user_args,
            'include_out': include_out,
            'randomization': randomization,
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
                'group_name': group_name,
                'out_group_name': out_group_name,
                'min_abundance': min_abundance,
                'name': group_name
            }
        ]

        if include_out:
            params['run_cfgs'].append({
                'factor': factor,
                'group': out_group,
                'out_group': group,
                'group_name': out_group_name,
                'out_group_name': group_name,
                'min_abundance': min_abundance,
                'name': out_group_name
                }
            )

        if len(errors_list) > 0:
            send_error_as_email(params, '\n'.join(errors_list))
            # TODO: Return form with error marked
            self.redirect('/')
            return

        pipeline = RunPipeline(params, inputs)
        pipeline.start()

        self.redirect('/')


def validate_inputs(params, inputs):
    mapping_file = inputs['mapping_file']
    factor = params['factor']
    group = params['group']

    labels = mapping_file[0].strip().strip('#').split(run_config.DELIM)

    errors_list = list()

    if 'SampleID' not in labels:
        errors_list.append('"SampleID" not in the headers of the sample ' +
                           '<-> group info file')
    if factor not in labels:
        errors_list.append(('"%s" not in the headers of the sample <-> ' +
                            'group info file') % factor)
    mapping_dict = get_categ_samples_dict(mapping_file, factor)
    for l in group:
        if l not in mapping_dict.keys():
            errors_list.append(
                'Interest group label %s is not in groupfile' % l
            )
    else:
        out_group = mapping_dict.keys()
        [out_group.remove(l) for l in group]

    if params['randomization'] not in ['row', 'column', 'table']:
        errors_list.append('Randomization option not valid')
    return (errors_list, mapping_dict, out_group)
