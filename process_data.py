from google.appengine.runtime import DeadlineExceededError

from itertools import compress
from biom.parse import parse_biom_table

from process_results import ProcessResultsPipeline

import datetime
import collections
import logging
import numpy

from storage import Results
from email_results import send_error_as_email, send_results_as_email
from randomize_data import randomize
import run_config

import pipeline
import pipeline.common


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        ntimes = int(params['ntimes'])

        true_res = run_data(inputs['mapping_dict'],
                            parse_biom_table(inputs['data']),
                            params['run_cfgs'])
        processing = []
        if ntimes <= run_config.MAX_NUM_PARALLEL:
            pipes = [1 for i in range(ntimes)]
        else:
            # list of number of times to run in each pipe
            pipes = [ntimes/run_config.MAX_NUM_PARALLEL
                     for i in range(run_config.MAX_NUM_PARALLEL)]
            left_over = ntimes % run_config.MAX_NUM_PARALLEL
            for i in range(left_over):
                pipes[i] += 1
        logging.info('Starting %d parallel tasks for randomized data',
                     len(pipes))
        for i in range(len(pipes)):
            res = yield RunRandomDataPipeline(inputs, params, pipes[i], i + 1,
                                              len(pipes))
            processing.append(res)
        with pipeline.InOrder():
            yield pipeline.common.Ignore(*processing)
            yield ProcessResultsPipeline(params, inputs, true_res)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]

        if self.was_aborted:
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(params, error)
        else:
            attachments = self.outputs.default.value
            send_results_as_email(params, attachments)
            Results.delete_entries(self.root_pipeline_id)
        self.cleanup()


class RunRandomDataPipeline(pipeline.Pipeline):
    def run(self, inputs, params, num, process_id, num_processes):
        logging.info('Pipeline %d of %d starting run of %d randomizations',
                     process_id, num_processes, num)
        start = datetime.datetime.now()
        res = {cfg['name']: dict() for cfg in params['run_cfgs']}
        parsed_data = parse_biom_table(inputs['data'])
        for i in range(num):
            try:
                logging.info('Pipeline %d of %d running run %d of %d',
                             process_id, num_processes, i + 1, num)
                random_mapping, random_data = randomize(inputs['mapping_dict'],
                                                        parsed_data,
                                                        params['random_opt'])
                collate_result(res, run_data(random_mapping, random_data,
                                             params['run_cfgs']))
                now = datetime.datetime.now()
                # assume the next run will take the average of completed runs
                if ((now - start) * (i + 2)) // (i + 1) > \
                   run_config.MAX_RUNNING_TIME:
                    logging.info('Pipeline %d of %d starting child process ' +
                                 'to avoid deadline',
                                 process_id, num_processes)
                    Results.add_entry(self.root_pipeline_id,
                                      self.pipeline_id, res)
                    yield RunRandomDataPipeline(inputs, params, num - i - 1,
                                                process_id, num_processes)
                    break
            except DeadlineExceededError:
                logging.info('Pipeline %d of %d ran out of time; ' +
                             'starting child process',
                             process_id, num_processes)
                Results.add_entry(self.root_pipeline_id, self.pipeline_id, res)
                yield RunRandomDataPipeline(inputs, params, num - i - 1,
                                            process_id, num_processes)
                break
            if i == num - 1:
                logging.info('Pipeline %s processed all runs',
                             self.pipeline_id)
                Results.add_entry(self.root_pipeline_id, self.pipeline_id, res)
                return


def collate_result(compiled, result):
    for cfg in result:
        for threshold, otus in result[cfg].items():
            if threshold in compiled[cfg]:
                compiled[cfg][threshold].update(otus)
            else:
                compiled[cfg][threshold] = collections.Counter(otus)


def run_data(mapping, data, run_cfgs):
    res = dict()
    for cfg in run_cfgs:
        res[cfg['name']] = get_core(mapping, data, cfg['group'])
    return res


def get_core(mapping, data, group):
    # a table of presence/absence data for just the interest group samples
    interest = data.filterSamples(lambda values, id, md: id in mapping[group])
    interest_samples = len(interest.SampleIds)
    presence_counts = interest.transformSamples(
        lambda l, id, md: numpy.array([v > 0 for v in l])
    ).sum(axis='observation')
    presence_fracs = [float(count) / interest_samples
                      for count in presence_counts]
    otus = [observation[2]['taxonomy']
            for observation in interest.iterObservations()]
    return {int(frac * 100):
            list(compress(otus, [presence >= frac
                                 for presence in presence_fracs]))
            for frac in run_config.FRACS}
