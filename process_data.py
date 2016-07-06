from google.appengine.runtime import DeadlineExceededError

from itertools import compress
from biom.parse import parse_biom_table

from process_results import ProcessResultsPipeline

import random
import datetime
import collections
import logging

from storage import Results
from email_results import send_error_as_email, send_results_as_email

import pipeline
import pipeline.common

# How many parallel pipes to have processing the randomized data
MAX_NUM_PARALLEL = 8
# If it looks like the task will run longer than this to do another itteration
# start a new task
MAX_RUNNING_TIME = datetime.timedelta(minutes=9)


class RunPipeline(pipeline.Pipeline):
    def run(self, params, inputs):
        logging.info('Starting run')
        NTIMES = int(params['ntimes'])

        true_res = run_data(inputs['data'], inputs['mapping_dict'],
                            params['run_cfgs'])
        processing = []
        if NTIMES <= MAX_NUM_PARALLEL:
            pipes = [1 for i in range(NTIMES)]
        else:
            # list of number of times to run in each pipe
            pipes = [NTIMES/MAX_NUM_PARALLEL for i in range(MAX_NUM_PARALLEL)]
            left_over = NTIMES % MAX_NUM_PARALLEL
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
            yield ProcessResultsPipeline(params, true_res)

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
        for i in range(num):
            try:
                logging.info('Pipeline %d of %d running run %d of %d',
                             process_id, num_processes, i + 1, num)
                randomized_mapping = shuffle_dicts(inputs['mapping_dict'])
                collate_result(res, run_data(inputs['data'],
                                             randomized_mapping,
                                             params['run_cfgs']))
                now = datetime.datetime.now()
                # assume the next run will take the average of completed runs
                if ((now - start) * (i + 2)) // (i + 1) > MAX_RUNNING_TIME:
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


def run_data(data, mapping, run_cfgs):
    res = dict()
    rich_data = parse_biom_table(data)
    for cfg in run_cfgs:
        res[cfg['name']] = get_core(mapping, rich_data, cfg['group'])
    return res


FRACS = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]


def get_core(mapping, data, group):
    # a table of presence/absence data for just the interest group samples
    interest = data.filterSamples(lambda values, id, md: id in mapping[group])
    interest_samples = len(interest.SampleIds)
    presence_counts = interest.reduce(lambda s, v: s + (v > 0),
                                      axis='observation')
    presence_fracs = [float(count) / interest_samples
                      for count in presence_counts]
    otus = [observation[2]['taxonomy']
            for observation in interest.iterObservations()]
    return {int(frac * 100):
            list(compress(otus, [presence > frac
                                 for presence in presence_fracs]))
            for frac in FRACS}


def shuffle_dicts(d):           # takes dictionary
    # If items(), keys(), values() are called with no intervening
    # modifications to the dictionary, the lists will directly correspond.
    keys = d.keys()
    values, lengths = get_values_from_dict(d)
    random.shuffle(values)
    new_values = divide_list(values, lengths)
    return dict(zip(keys, new_values))


# this can be changed later to allow splitting in more than two groups
def divide_list(a, lengths):
    # a[start:end] # items start through end-1
    return a[:lengths[0]], a[lengths[0]:]


def get_values_from_dict(a):
    values = list()
    lengths = list()
    for v in a.values():
        values.extend(v)
        lengths.append(len(v))  # sizes of original lists
    return values, lengths
