from compute_core_microbiome import exec_core_microb_cmd
from process_results import ProcessResultsPipeline

import random

from storage import Result_RandomDict, clean_storage
from email_results import send_error_as_email, send_results_as_email

import pipeline
import pipeline.common

# How many parallel pipes to have processing the randomized data
MAX_NUM_PARALLEL = 50


class RunPipeline(pipeline.Pipeline):
    def run(self, params):
        data = params['data']
        mapping_file = params['mapping_file']
        mapping_dict = params['mapping_dict']
        factor = params['factor']
        group = params['group']
        out_group = params['out_group']
        DELIM = params['delim']
        NTIMES = int(params['ntimes'])

        core = yield RunDataPipeline(data, mapping_file, factor, group, DELIM)
        out = yield RunDataPipeline(data, mapping_file, factor, out_group,
                                    DELIM)
        processing = []
        if NTIMES <= MAX_NUM_PARALLEL:
            pipes = [1 for i in range(NTIMES)]
        else:
            # list of number of times to run in each pipe
            pipes = [NTIMES/MAX_NUM_PARALLEL for i in range(MAX_NUM_PARALLEL)]
            left_over = NTIMES % MAX_NUM_PARALLEL
            for i in range(left_over):
                pipes[i] += 1
        for num_in_task in pipes:
            res = yield RunRandomDataPipeline(data, mapping_dict,
                                              factor, group, out_group, DELIM,
                                              num_in_task)
            processing.append(res)
        with pipeline.InOrder():
            yield pipeline.common.Ignore(*processing)
            yield ProcessResultsPipeline(params, core, out)

    def finalized(self):
        params = self.args[0]
        timestamp = params['timestamp']
        user_args = params['user_args']
        to_email = params['to_email']

        if self.was_aborted:
            error = 'An unknown error has occured. Please try again. ' +\
                    'If this occurs again please contact the developers'
            send_error_as_email(timestamp, user_args, error, to_email)
        else:
            results_string = self.outputs.default.value[0]
            tree = self.outputs.default.value[1]
            send_results_as_email(timestamp, user_args, results_string, tree,
                                  to_email)
        clean_storage(self.root_pipeline_id)


class RunRandomDataPipeline(pipeline.Pipeline):
    def run(self, data, mapping_dict, factor, group, out_group, DELIM,
            num):
        processing = []
        # To prevent spawning too many tasks at one time
        with pipeline.InOrder():
            for i in range(num):
                randomized_mapping = convert_shuffled_dict_to_str(
                    shuffle_dicts(mapping_dict), factor)
                core = yield RunDataPipeline(data, randomized_mapping,
                                             factor, group, DELIM)
                out = yield RunDataPipeline(data, randomized_mapping,
                                            factor, out_group, DELIM,
                                            out_group=True)
                write = yield WriteRandomResultPipeline(self.pipeline_id,
                                                        core, out)
                processing.append(write)
            # Wait till everything is done
        yield pipeline.common.Ignore(*processing)


class RunDataPipeline(pipeline.Pipeline):
    def run(self, data, mapping, factor, group, DELIM, out_group=False):
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


class WriteRandomResultPipeline(pipeline.Pipeline):
    def run(self, run_id, core, out):
        return Result_RandomDict.add_entry(self.root_pipeline_id,
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
