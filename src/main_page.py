import webapp2
from time import localtime, strftime
import cPickle

import run_config
from process_data import RunPipeline
from parse_inputs import parse_inputs


class MainPage(webapp2.RequestHandler):
    def get(self):
        template = run_config.JINJA_ENVIRONMENT.get_template('index.html')
        self.response.write(template.render())

    def post(self):
        # Used to measure total processing time
        timestamp = strftime('%a-%d-%b-%Y-%I:%M:%S-%p', localtime())

        run_name = self.request.get('name')
        data = self.request.get('datafile')
        mapping_file = self.request.get('groupfile').split('\n')

        factor = self.request.get('factor')
        group = map(lambda s: s.strip(), self.request.get('group').split(','))
        group_name = self.request.get('group_name')
        out_group_name = self.request.get('out_group_name')

        p_val_adj = self.request.get('pvaladjmethod')
        include_out = bool(self.request.get('include_out'))
        max_p = float(self.request.get('max_p'))
        min_abundance = float(self.request.get('min_abundance'))

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\nMax Pval: %s\n\n\n')
                     % (factor, ', '.join(group), p_val_adj, max_p))

        params = {
            'run_name': run_name,
            'factor': factor,
            'group': group,
            'to_email': to_email,
            'timestamp': timestamp,
            'user_args': user_args,
            'include_out': include_out,
            'max_p': max_p,
            'min_abundance': min_abundance,
        }

        errors_list, mapping_dict, out_group, filtered_data \
            = parse_inputs(params, mapping_file, data)

        params['out_group'] = out_group
        inputs = {
            'mapping_dict': mapping_dict,
            # Must be packed up to be able to pass to the process_data task
            'filtered_data': cPickle.dumps(filtered_data),
        }

        params['run_cfgs'] = [
            {
                'factor': factor,
                'group': group,
                'out_group': out_group,
                'group_name': group_name,
                'out_group_name': out_group_name,
                'min_abundance': min_abundance,
                'name': group_name,
                'max_p': max_p,
                'p_val_adj': p_val_adj,
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
                'name': out_group_name,
                'max_p': max_p,
                'p_val_adj': p_val_adj,
                }
            )

        template = run_config.JINJA_ENVIRONMENT.get_template('result.html')
        if len(errors_list) > 0:
            self.response.write(template.render(errors=errors_list,
                                                sucess=False))
            return

        pipeline = RunPipeline(params, inputs)
        pipeline.start()

        self.response.write(template.render(user_email=to_email, sucess=True))
