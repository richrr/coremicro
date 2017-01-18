# Copyright 2016, 2017 Richard Rodrigues, Nyle Rodgers, Mark Williams,
# Virginia Tech
#
# This file is part of Coremic.
#
# Coremic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Coremic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Coremic. If not, see <http://www.gnu.org/licenses/>.


import webapp2
from time import localtime, strftime
import cPickle

import web_config
from run_pipeline import RunPipeline
from ..core.parse_inputs import parse_inputs


class MainPage(webapp2.RequestHandler):
    def get(self):
        template = web_config.JINJA_ENVIRONMENT.get_template('index.html')
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
        min_frac = float(self.request.get('min_frac'))
        min_abundance = float(self.request.get('min_abundance'))

        to_email = self.request.get('email')

        user_args = (('You selected the following parameters:' +
                      '\nFactor: %s\nGroup: %s\n' +
                      'Pval correction: %s\nMax Pval: %s\n\n\n')
                     % (factor, ', '.join(group), p_val_adj, max_p))

        params = {
            'run_name': run_name,
            'to_email': to_email,
            'timestamp': timestamp,
            'user_args': user_args,
            'factor': factor,
            'group': group,
            'max_p': max_p,
        }

        errors_list, mapping_dict, out_group, filtered_data, original_otus \
            = parse_inputs(params, mapping_file, data)

        params['user_args'] += (
            'OTUs before combining duplicates: %s\n' +
            'OTUs after combining duplicates: %s\n\n\n') % (
                original_otus,
                len(filtered_data.ObservationIds)
            )

        params['out_group'] = out_group
        inputs = {
            'mapping_dict': mapping_dict,
            # Must be packed up to be able to pass to the process_data task
            'filtered_data': cPickle.dumps(filtered_data),
        }

        params['run_cfgs'] = [
            {
                'run_name': run_name,
                'factor': factor,
                'group': group,
                'out_group': out_group,
                'group_name': group_name,
                'out_group_name': out_group_name,
                'min_abundance': min_abundance,
                'name': group_name,
                'max_p': max_p,
                'min_frac': min_frac,
                'p_val_adj': p_val_adj,
            }
        ]

        if include_out:
            params['run_cfgs'].append({
                'run_name': run_name,
                'factor': factor,
                'group': out_group,
                'out_group': group,
                'group_name': out_group_name,
                'out_group_name': group_name,
                'min_abundance': min_abundance,
                'name': out_group_name,
                'max_p': max_p,
                'min_frac': min_frac,
                'p_val_adj': p_val_adj,
                }
            )

        template = web_config.JINJA_ENVIRONMENT.get_template('result.html')
        if len(errors_list) > 0:
            self.response.write(template.render(errors=errors_list,
                                                sucess=False))
            return

        pipeline = RunPipeline(params, inputs)
        pipeline.start()

        self.response.write(template.render(user_email=to_email, sucess=True))
