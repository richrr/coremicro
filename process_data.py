from compute_core_microbiome import exec_core_microb_cmd
from process_results import ProcessResultsPipeline

import random

from storage import OriginalBiom, Result_RandomDict

import pipeline
import pipeline.common

MIN_TASKS_PER_PIPE = 50


class RunPipeline(pipeline.Pipeline):
    def run(self, key):
        user_args, to_email, p_val_adj, DELIM, NTIMES, data, \
            mapping, factor, group, out_group, OUTPFILE, \
            mapping_dict = OriginalBiom.get_params(key)

        core = yield RunDataPipeline(key, data, mapping, factor, group, DELIM)
        out = yield RunDataPipeline(key, data, mapping, factor, out_group,
                                    DELIM)
        results = []
        for i in range(int(int(NTIMES)/MIN_TASKS_PER_PIPE)):
            # TODO: Don't round down to nearest 50!
            res = yield RunRandomDataPipeline(key, data, mapping_dict,
                                              factor, group, out_group, DELIM,
                                              MIN_TASKS_PER_PIPE)
            results.append(res)
        with pipeline.InOrder():
            yield pipeline.common.Ignore(*results)
            yield ProcessResultsPipeline(key, core, out)


class RunRandomDataPipeline(pipeline.Pipeline):
    def run(self, key, data, mapping_dict, factor, group, out_group, DELIM,
            num):
        processing = []
        for i in range(num):
            randomized_mapping = convert_shuffled_dict_to_str(
                shuffle_dicts(mapping_dict), factor)
            core = yield RunDataPipeline(key, data, randomized_mapping,
                                         factor, group, DELIM)
            core_write = yield WriteRandomResults(key, core,
                                                  out_group=False)
            processing.append(core_write)
            out = yield RunDataPipeline(key, data, randomized_mapping,
                                        factor, out_group, DELIM)
            out_write = yield WriteRandomResults(key, out, out_group=True)
            processing.append(out_write)
        # Wait till everything is done
        yield pipeline.common.Ignore(*processing)


class WriteRandomResults(pipeline.Pipeline):
    def run(self, key, results, out_group=False):
        return Result_RandomDict.add_entry(key, results, out_group)


class RunDataPipeline(pipeline.Pipeline):
    def run(self, key, data, mapping, factor, group, DELIM):
        result = exec_core_microb_cmd(data, 'dir', mapping, factor, group)
        # compile original results. The key is frac_threshold,
        # value is a list of unique otus
        result_frac_thresh_otus_dict = dict()
        # return the items in sorted order
        for frac_thresh, core_OTUs_biom in (
                sorted(result['frac_thresh_core_OTUs_biom'].items(),
                       key=lambda (key, value): int(key))):
            OTUs, biom = core_OTUs_biom  # this is a tuple of otus and biom
            result_frac_thresh_otus_dict[frac_thresh] = compile_results(OTUs,
                                                                        DELIM)

        return result_frac_thresh_otus_dict


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
