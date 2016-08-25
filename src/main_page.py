import webapp2
import jinja2
import os
from time import localtime, strftime

from process_data import RunPipeline
from parse_inputs import get_categ_samples_dict, read_table
import run_config


JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__),
                                                '../static/')),
    extensions=['jinja2.ext.autoescape'],
    autoescape=True)


class MainPage(webapp2.RequestHandler):
    def get(self):
        template = JINJA_ENVIRONMENT.get_template('index.html')
        self.response.write(template.render())

    def post(self):
        name = self.request.get('name')
        data = self.request.get('datafile')
        mapping_file = self.request.get('groupfile').split('\n')

        factor = self.request.get('factor')
        group = map(lambda s: s.strip(), self.request.get('group').split(','))
        group_name = self.request.get('group_name')
        out_group_name = self.request.get('out_group_name')
        timestamp = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())

        p_val_adj = self.request.get('pvaladjmethod')
        include_out = bool(self.request.get('include_out'))
        min_abundance = float(self.request.get('min_abundance'))

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\n\n\n')
                     % (factor, group, p_val_adj))

        params = {
            'name': name,
            'factor': factor,
            'group': group,
            'p_val_adj': p_val_adj,
            'to_email': to_email,
            'timestamp': timestamp,
            'user_args': user_args,
            'include_out': include_out,
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
            template = JINJA_ENVIRONMENT.get_template('result.html.jinja')
            self.response.write(template.render(errors=errors_list,
                                                sucess=False))
            return

        pipeline = RunPipeline(params, inputs)
        pipeline.start()

        template = JINJA_ENVIRONMENT.get_template('result.html.jinja')
        self.response.write(template.render(user_email=to_email, sucess=True))


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

    try:
        out_group = mapping_dict.keys()
        [out_group.remove(l) for l in group]
    except ValueError:
        pass                    # already handled previously

    try:
        read_table(inputs['data'])
    except ValueError as e:
        errors_list.append('Datafile could not be read: %s' % e.message)
    return (errors_list, mapping_dict, out_group)
