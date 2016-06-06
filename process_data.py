import webapp2

from compute_core_microbiome import exec_core_microb_cmd

import random

from storage import Result_RandomDict, Result_TrueDict, OriginalBiom


class ProcessData(webapp2.RequestHandler):
    def post(self):
        key = self.request.get('otu_table_biom_key')
        mode = self.request.get('mode')

        user_args, to_email, p_val_adj, DELIM, NTIMES, otu_table_biom, \
            mapping_info_list, factor, group, out_group, OUTPFILE, \
            categ_samples_dict = OriginalBiom.get_params(key)

        if mode == 'true':
            res = run_data(otu_table_biom, OUTPFILE,
                           mapping_info_list, factor, group, DELIM)
            Result_TrueDict.add_entry(key, res)
        elif mode == 'out':
            res = run_data(otu_table_biom, OUTPFILE,
                           mapping_info_list, factor, out_group, DELIM)
            Result_TrueDict.add_entry(key, res, out_group=True)
        elif mode == 'random':
            random_info_lists = randomize_info(factor, 50, categ_samples_dict)

            for rand_list in random_info_lists:
                Result_RandomDict.add_entry(key,
                                            run_data(otu_table_biom, OUTPFILE,
                                                     rand_list, factor, group,
                                                     DELIM))
                Result_RandomDict.add_entry(key,
                                            run_data(otu_table_biom, OUTPFILE,
                                                     rand_list, factor,
                                                     out_group, DELIM),
                                            out_group=True)


def run_data(otu_table_biom, o_dir, mapping_info_list, c, group, DELIM):
    result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list,
                                  c, group)
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


def randomize_info(factor, NTIMES, categ_samples_dict):
    print "Creating shuffled dicts"
    random_info_lists = []
    for i in range(NTIMES):
        shuffled_dict = shuffle_dicts(categ_samples_dict)
        rand_mapping_info_list = convert_shuffled_dict_to_str(shuffled_dict,
                                                              factor)
        random_info_lists.append(rand_mapping_info_list)
    return random_info_lists


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
