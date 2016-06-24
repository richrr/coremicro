from google.appengine.runtime import DeadlineExceededError

from compute_core_microbiome import exec_core_microb_cmd
from process_results import ProcessResultsPipeline

import random
import datetime
import collections
import logging

from storage import Result_RandomDict, clean_storage
from email_results import send_error_as_email, send_results_as_email

import pipeline
import pipeline.common

# How many parallel pipes to have processing the randomized data
MAX_NUM_PARALLEL = 20
# If it looks like the task will run longer than this to do another itteration
# start a new task
MAX_RUNNING_TIME = datetime.timedelta(minutes=9)


class RunPipeline(pipeline.Pipeline):
    def run(self, params):
        logging.info('Starting run')
        data = params['data']
        mapping_file = params['mapping_file']
        mapping_dict = params['mapping_dict']
        factor = params['factor']
        group = params['group']
        out_group = params['out_group']
        DELIM = params['delim']
        NTIMES = int(params['ntimes'])

        core = run_data(data, mapping_file, factor, group, DELIM)
        out = run_data(data, mapping_file, factor, out_group, DELIM)
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
        for num_in_task in pipes:
            res = yield RunRandomDataPipeline(data, mapping_dict,
                                              factor, group, out_group, DELIM,
                                              num_in_task)
            processing.append(res)
        with pipeline.InOrder():
            yield pipeline.common.Ignore(*processing)
            yield ProcessResultsPipeline(params, core, out)

    def finalized(self):
        logging.info('Finalizing task')
        params = self.args[0]
        timestamp = params['timestamp']
        user_args = params['user_args']
        to_email = params['to_email']
        name = params['name']

        if self.was_aborted:
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(timestamp, user_args, error, name, to_email)
        else:
            results_string = self.outputs.default.value[0]
            tree = self.outputs.default.value[1]
            send_results_as_email(timestamp, user_args, results_string, tree,
                                  name, to_email)
        clean_storage(self.root_pipeline_id)
        self.cleanup()


class RunRandomDataPipeline(pipeline.Pipeline):
    def run(self, data, mapping_dict, factor, group, out_group, DELIM, num):
        logging.info('Pipeline %s starting run of %d randomizations',
                     self.pipeline_id, num)
        start = datetime.datetime.now()
        core_comp = dict()
        out_comp = dict()
        for i in range(num):
            try:
                logging.info('Pipeline %s running run %d of %d',
                             self.pipeline_id, i + 1, num)
                randomized_mapping = convert_shuffled_dict_to_str(
                    shuffle_dicts(mapping_dict), factor)
                add_result(core_comp, run_data(data, randomized_mapping,
                                               factor, group, DELIM))
                add_result(out_comp, run_data(data, randomized_mapping,
                                              factor, out_group, DELIM))
                now = datetime.datetime.now()
                # assume the next run will take the average of completed runs
                if ((now - start) * (i + 2)) // (i + 1) > MAX_RUNNING_TIME:
                    logging.info(
                        'Pipeline %s starting child process to avoid deadline',
                        self.pipeline_id)
                    write_random_result(self.root_pipeline_id,
                                        self.pipeline_id, core_comp, out_comp)
                    yield RunRandomDataPipeline(data, mapping_dict,
                                                factor, group,
                                                out_group, DELIM,
                                                num - i - 1)
                    break
            except DeadlineExceededError:
                logging.info(
                    'Pipeline %s ran out of time; starting child process',
                    self.pipeline_id)
                write_random_result(self.root_pipeline_id, self.pipeline_id,
                                    core_comp, out_comp)
                yield RunRandomDataPipeline(data, mapping_dict, factor,
                                            group, out_group, DELIM,
                                            num - i - 1)
                break
            if i == num - 1:
                logging.info('Pipeline %s processed all runs',
                             self.pipeline_id)
                write_random_result(self.root_pipeline_id, self.pipeline_id,
                                    core_comp, out_comp)
                return


def add_result(compiled, result):
    for threshold, otus in result.items():
        if threshold in compiled:
            compiled[threshold].update(otus)
        else:
            compiled[threshold] = collections.Counter(otus)


def run_data(data, mapping, factor, group, DELIM):
    result = exec_core_microb_cmd(data, 'dir', mapping, factor, group)
    # compile original results. The key is frac_threshold,
    # value is a list of unique otus
    compiled = dict()
    # return the items in sorted order
    for frac_thresh, core_OTUs_biom in (
            sorted(result['frac_thresh_core_OTUs_biom'].items(),
                   key=lambda (key, value): int(key))):
        OTUs, biom = core_OTUs_biom  # this is a tuple of otus and biom
        compiled[frac_thresh] = compile_results(OTUs, DELIM)

    return compiled


def write_random_result(root_id, run_id, core, out):
    logging.info('Writing results from pipeline %s', run_id)
    return Result_RandomDict.add_entry(root_id,
                                       run_id, core, out)


# get unique elements from last column (otus)
def compile_results(otus, DELIM):
    taxon_list = list()         # this may be a json or list
    result_ = otus              # json.loads(otus)
    for l in result_:
        l = l.strip()
        contents = l.split(DELIM)
        if '#' in l or not l:
            continue
        taxon_list.append(contents[1])
    return list(set(taxon_list))


def convert_shuffled_dict_to_str(DICT, categ):
    file_str_list = list()
    file_str = "%s\t%s\n" % ('#SampleID', categ)
    file_str_list.append(file_str)

    for k, v in DICT.items():
        for i in v:
            f_str = "%s\t%s\n" % (i, k)
            file_str_list.append(f_str)
    return file_str_list


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
