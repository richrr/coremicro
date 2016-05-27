import webapp2

from google.appengine.ext import ndb

from compute_core_microbiome import exec_core_microb_cmd

import sys
import random
from utils import list_to_dict

from storage import Result_RandomDict, Result_TrueDict
from utilities import compile_results, get_required_params_from_orig_dict
from email_results import send_results_as_email


class ProcessData(webapp2.RequestHandler):
    def post(self):
        otu_table_biom_o = self.request.get('otu_table_biom_key')
        mode = self.request.get('mode')

        # get rid of this orig biom and use only the method
        # also move the true results to a separate task to avoid redoing work
        user_args, to_email, p_val_adj, DELIM, NTIMES, otu_table_biom, \
            g_info_list, factor, group, out_group, OUTPFILE = (
                get_required_params_from_orig_dict(otu_table_biom_o))

        indx_sampleid, indx_categ, errors_list = validate_inputs(
            'ndb_custom_key', 'user_args', otu_table_biom, factor, group,
            out_group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE,
            to_email)
        # temp hack since blobstore randomly swaps file order during upload
        if len(errors_list) > 0:
            tmp = g_info_list
            g_info_list = otu_table_biom
            otu_table_biom = tmp
            # retry with swapped files
            indx_sampleid, indx_categ, errors_list = validate_inputs(
                'ndb_custom_key', 'user_args', otu_table_biom, factor, group,
                out_group, g_info_list, p_val_adj, DELIM, int(NTIMES),
                OUTPFILE, to_email)
            if len(errors_list) > 0:  # just give up on this
                print '\n'.join(errors_list)
                # put code here so that the code doesn't run further
                sys.exit(0)
            else:
                print 'Swapping files worked!'
        else:
            print 'No file swapping needed!'

        # this is to query all entries in this run
        ndb_custom_key = OUTPFILE + '~~~~' + factor + '~~~~' + group
        ndb_custom_key_out = OUTPFILE + '~~~~' + factor + '~~~~' + out_group

        if mode == 'true':
            run_true_data(OUTPFILE, otu_table_biom, g_info_list, factor,
                          group, DELIM, ndb_custom_key)
        elif mode == 'out':
            run_true_data(OUTPFILE, otu_table_biom, g_info_list, factor,
                          group, DELIM, ndb_custom_key_out)
        elif mode == 'random':
            calc_significance(indx_sampleid, indx_categ, errors_list,
                              otu_table_biom, factor, group, g_info_list,
                              p_val_adj, DELIM, int(NTIMES), OUTPFILE,
                              to_email)


# run core microbiome on original data
def run_true_data(OUTPFILE, otu_table_biom, mapping_info_list, c, group, DELIM,
                  ndb_custom_key):
    o_dir = 'true_result' + OUTPFILE
    result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list,
                                  c, group)

    # compile original results. The key is frac_threshold,
    # value is a list of unique otus
    true_result_frac_thresh_otus_dict = dict()
    # return the items in sorted order
    for frac_thresh, core_OTUs_biom in (
            sorted(result['frac_thresh_core_OTUs_biom'].items(),
                   key=lambda (key, value): int(key))):
        OTUs, biom = core_OTUs_biom  # this is a tuple of otus and biom
        true_result_frac_thresh_otus_dict[frac_thresh] = compile_results(OTUs,
                                                                         DELIM)

    fatherrestrue_idx = 'fatherresultstrue' + ndb_custom_key
    Result_TrueDict(parent=ndb.Key(Result_TrueDict, fatherrestrue_idx),
                    idx=ndb_custom_key,
                    true_results=true_result_frac_thresh_otus_dict).put()
    print ("Processed %s fraction thresholds for true data" %
           str(len(true_result_frac_thresh_otus_dict)))


def calc_significance(indx_sampleid, indx_categ, errors_list,
                      otu_table_biom, c, group, mapping_info_list,
                      p_val_adj, DELIM, NTIMES, OUTPFILE, to_email):

    # global ndb_custom_key
    # this is to query all entries in this run
    ndb_custom_key = OUTPFILE + '~~~~' + c + '~~~~' + group
    user_args = (('You selected the following parameters:\nFactor: %s\n' +
                  'Group: %s\nPval correction: %s\n' +
                  '# of randomizations: %s\n\n\n')
                 % (c, group, p_val_adj, NTIMES))

    categ_samples_dict = list_to_dict(mapping_info_list, DELIM, ',',
                                      'current', indx_categ, indx_sampleid)

    if not check_map_file_has_two_groups(categ_samples_dict):
        errors_list.append('\nERROR: Following code divides samples in ' +
                           '>TWO groups. Change the mapping file to only '
                           'have two groups (e.g. A vs D)\n')

    # email the error and quit, no point to continue further
    if len(errors_list) > 0:
        send_results_as_email(ndb_custom_key, user_args,
                              '\n'.join(errors_list), to_email)
        # put code here so that the code doesn't run further

    print "Creating shuffled dicts"
    for i in range(NTIMES):
        shuffled_dict = shuffle_dicts(categ_samples_dict)
        ndb_custom_key_entry = ndb_custom_key + '~~~~' + str(i)
        shuffle_dict_coremic_serial_dict_datastore(shuffled_dict,
                                                   ndb_custom_key,
                                                   ndb_custom_key_entry,
                                                   otu_table_biom)


def shuffle_dict_coremic_serial_dict_datastore(shuffled_dict,
                                               ndb_custom_key,
                                               ndb_custom_key_entry,
                                               otu_table_biom):
    '''
    get the randomized dict and run core microbiome
    '''

    local_dict_frac_thresh_otus = dict()
    # Zen-outputMon-07-Mar-2016-01:36:13-AM~~~~Plant~~~~Sw~~~~2
    OUTPFILE, c, group, rand_iter_numb = ndb_custom_key_entry.split('~~~~')
    rand_mapping_info_list = convert_shuffled_dict_to_str(
        shuffled_dict, c)
    rand_o_dir = rand_iter_numb + OUTPFILE
    result = exec_core_microb_cmd(otu_table_biom, rand_o_dir,
                                  rand_mapping_info_list, c, group)

    # arrange the results to look pretty
    # return the items in sorted order
    for r_frac_thresh, r_core_OTUs_biom in (
            sorted(result['frac_thresh_core_OTUs_biom'].items(),
                   key=lambda (key, value): int(key))):
        r_OTUs, r_biom = r_core_OTUs_biom

        ndb_custom_key_r_frac_thres = (ndb_custom_key + '~~~~' +
                                       r_frac_thresh)

        if ndb_custom_key_r_frac_thres in local_dict_frac_thresh_otus:
            print "Why do you have same fraction thresholds repeating?"
        else:
            local_dict_frac_thresh_otus[ndb_custom_key_r_frac_thres] = (
                r_OTUs)

    fatherres_idx = 'fatherresults' + ndb_custom_key
    Result_RandomDict(parent=ndb.Key(Result_RandomDict, fatherres_idx),
                      idx=ndb_custom_key,
                      otus=local_dict_frac_thresh_otus).put()


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
    for i in a.values():
        v = i.split(',')
        values.extend(v)
        lengths.append(len(v))  # sizes of original lists
    return values, lengths


# takes iternumb which is number of random iteration and dictionary
def check_map_file_has_two_groups(a):
    values, lengths = get_values_from_dict(a)
    if len(lengths) > 2:
        return False
    return True


def validate_inputs(ndb_custom_key, user_args, otu_table_biom, c, group,
                    out_group, mapping_info_list, p_val_adj, DELIM, NTIMES,
                    OUTPFILE, to_email):
    # find index of SampleID and category to be summarized. e.g. swg or non-swg
    labels = mapping_info_list[0].split(DELIM)
    indx_sampleid = indx_categ = ''

    errors_list = list()

    try:  # this is from the mapping file
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        errors_list.append("' not in the headers of the sample <-> " +
                           "group info file")
    try:
        indx_categ = labels.index(c)
    except ValueError:
        errors_list.append("'%s' not in the headers of the sample <-> " +
                           "group info file" % c)
    if int(NTIMES) % 50 != 0:
        errors_list.append("Number of randomizations requested is not " +
                           "multiple of 50. Kindly rerun job")
    return (indx_sampleid, indx_categ, errors_list)
