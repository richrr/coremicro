# Copyright 2016 Richard Rodrigues, Nyle Rodgers, Mark Williams, Virginia Tech
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
import logging
import pipeline
import cPickle
import base64
from datetime import datetime
from time import strptime, mktime

from send_email import send_email
from generate_graph import generate_graph
from core.process_data import process, format_results


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        # Unpack the data
        inputs['filtered_data'] = cPickle.loads(str(inputs['filtered_data']))

        attachments = list()
        for cfg in params['run_cfgs']:
            core = process(inputs, cfg)
            attachments += ([{
                'Content-Type': 'text/plain',
                'Filename': '%s_results_%s.tsv' % (cfg['name'],
                                                   cfg['run_name']),
                'content': base64.b64encode(format_results(core, cfg))}] +
                            generate_graph(inputs, cfg, core))
        elapsed_time = datetime.now() - datetime.fromtimestamp(mktime(
            strptime(params['timestamp'], '%a-%d-%b-%Y-%I:%M:%S-%p')))
        send_email(('Your data with name %s has been processed' %
                    params['run_name']),
                   '''Dear User:

Your data has been processed and is attached. Thanks for using this tool.

Please email us if you have any questions.

The Core Microbiome Team

%s
Elapsed time: %d.%06d seconds''' % (params['user_args'],
                                    elapsed_time.seconds,
                                    elapsed_time.microseconds),
                   params['to_email'], attachments)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]

        if self.was_aborted:    # There's probably a bug in the code
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            logging.warn(error)
            elapsed_time = datetime.now() - datetime.fromtimestamp(mktime(
                strptime(params['timestamp'], '%a-%d-%b-%Y-%I:%M:%S-%p')))
            send_email('There was an error in processing your data with ' +
                       'name %s' % params['run_name'],
                       '''Dear User:

There was an error in processing your data. The error is listed below.

Please email us if you have any questions.

The Core Microbiome Team

%s
Elapsed time %d.%06d seconds

%s''' % (params['user_args'], elapsed_time.seconds, elapsed_time.microseconds,
         error),
                       params['to_email'])
        self.cleanup()
