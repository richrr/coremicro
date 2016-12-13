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

import run_config
from send_email import send_email
from parse_inputs import get_data_summary
from probability import row_randomize_probability
from generate_graph import generate_graph
from correct_p_values import correct_pvalues


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        # Unpack the data
        inputs['filtered_data'] = cPickle.loads(str(inputs['filtered_data']))

        attachments = list()
        for cfg in params['run_cfgs']:
            results = get_signif_otus(inputs, cfg)
            attachments += format_results(results, params, cfg)
            attachments += generate_graph(params, inputs, cfg, results)
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


def get_signif_otus(inputs, cfg):
    """Takes the set of OTUs at each fractional threshold that satisfy the
    inclusion requirement that they be present at at least the threshold in the
    interest group and at less than that threshold across all samples,
    calculates the p-value for each OTU, and gives the list of OTUs found to be
    significant after correction for multiple testing"""
    data = get_data_summary(inputs['filtered_data'],
                            inputs['mapping_dict'],
                            cfg['group'],
                            cfg['min_abundance'])
    results = dict()
    core = {frac: [otu for otu in data.keys()
                   if (data[otu].i_present /
                       float(data[otu].n_interest) >= frac)]
            for frac in run_config.FRACS}
    for frac in core:
        logging.info('Starting frac %f' % frac)
        pvals = [row_randomize_probability(data[otu])
                 for otu in core[frac]]
        pvals_corrected = correct_pvalues(pvals, cfg['p_val_adj'])
        results[frac] = [{'otu': otu,
                          'pval': pvals[i],
                          'corrected_pval': pvals_corrected[i],
                          'threshold': frac}
                         for i, otu in enumerate(core[frac])
                         if pvals_corrected[i] <= cfg['max_p']]
    return results


def format_results(res, params, cfg):
    """Format the result data as a tsv
    """
    sign_results = (('OTU\tpval\t%s ' +
                     'corrected pval\tthreshold\n')
                    % cfg['p_val_adj'])
    for frac in reversed(sorted(res.keys())):
        for otu in res[frac]:
            sign_results += '%s\t%s\t%s\t%s\n' % (
                otu['otu'], otu['pval'],
                otu['corrected_pval'], int(frac * 100))
    if not run_config.IS_PRODUCTION:
        print "Results for configuration: " + cfg['name']
        print sign_results
    return [{
        'Content-Type': 'text/plain',
        'Filename': '%s_results_%s.tsv' % (cfg['name'], params['run_name']),
        'content': base64.b64encode(sign_results)
    }]
