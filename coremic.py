#!/usr/bin/env python
#
# Copyright 2007 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import urllib
import jinja2

import cgi
import datetime
import webapp2

import time

from google.appengine.ext import ndb
from google.appengine.api import users

from utils import *

import sys

import operator
from time import localtime, strftime

#import argparse
#import brewer2mpl

import math

import numpy
#import matplotlib.pyplot as plt
#import pylab
#from scipy import stats
#import matplotlib as mpl

import string
import random
#from subprocess import Popen, PIPE
import collections
#from multiprocessing import Pool, Lock
import itertools

from convert_biom import *
from compute_core_microbiome import *


from google.appengine.ext import blobstore
#from google.appengine.ext import db

from google.appengine.ext.webapp import blobstore_handlers

'''
from google.appengine.api import app_identity
'''
from mapreduce import base_handler
from mapreduce import mapreduce_pipeline
from mapreduce import operation as op
from mapreduce import shuffler
from mapreduce import context

import json

from pipeline import pipeline

from google.appengine.api import mail

import logging
from google.appengine.runtime import DeadlineExceededError
from google.appengine.runtime import apiproxy_errors
from google.appengine.api import taskqueue
from google.appengine.api import background_thread

from google.appengine.api.taskqueue import TaskRetryOptions

from google.appengine.ext.blobstore import BlobKey

#from google.appengine.ext import deferred

# to do the parallelism: async, futures or mapreduce
# https://cloud.google.com/appengine/docs/python/datastore/async#Working_with_the_Async_Datastore_API
# futures is for python 3, so cannot be used for my 2.7
# mapreduce should work; if it doesn't, try the async datastore api


#http://stackoverflow.com/questions/11849456/how-to-filter-datastore-data-before-mapping-to-cloud-storage-using-the-mapreduce
#http://stackoverflow.com/questions/23508116/appengine-mapreduce-how-to-filter-structuredproperty-while-using-datastore-input


## for every random dict entry, it has different thresholds
# the actual key is automatically generated
class Result_RandomDict(ndb.Model):
  idx = ndb.StringProperty() # the run id
  #frac_thresh = ndb.StringProperty() # the frac threshold along with the run id
  #entry_id = ndb.StringProperty() # the frac threshold and entry along with the run id
  otus = ndb.JsonProperty() # the core otus from the (shuffled) dictionary as a json
  #true_results = ndb.JsonProperty() # the core biom from the true dictionary as a json


class Result_TrueDict(ndb.Model):
  idx = ndb.StringProperty() # the run id
  #frac_thresh = ndb.StringProperty() # the frac threshold along with the run id
  #entry_id = ndb.StringProperty() # the frac threshold and entry along with the run id
  #otus = ndb.JsonProperty() # the core otus from the (shuffled) dictionary as a json
  true_results = ndb.JsonProperty() # the core biom from the true dictionary as a json


# the actual key is automatically generated
class OriginalBiom(ndb.Model):
  idx = ndb.StringProperty() # the run id
  biom = ndb.JsonProperty() # the core biom from the original input as a json
  params_str = ndb.JsonProperty()


# this cane be changed later to allow splitting in more than two groups
def divide_list(a, lengths):
    return a[:lengths[0]] , a[lengths[0]:]  # a[start:end] # items start through end-1





def shuffle_dict_coremic_serial_dict(entty, ndb_custom_key, otu_table_biom):
    '''
     get the randomized dict and run core microbiome
    '''    
    r_out_str = ''
    local_dict_frac_thresh_otus = dict()
    if entty['idx'] == ndb_custom_key:
        rand_mapping_info_dict = entty['dict']
        OUTPFILE, c, group, rand_iter_numb = entty['entry_id'].split('~~~~') #Zen-outputMon-07-Mar-2016-01:36:13-AM~~~~Plant~~~~Sw~~~~2
        rand_mapping_info_list = convert_shuffled_dict_to_str(rand_mapping_info_dict, c)
        rand_o_dir = rand_iter_numb + OUTPFILE
        result = exec_core_microb_cmd(otu_table_biom, rand_o_dir, rand_mapping_info_list, c, group)

        '''
         arrange the results to look pretty
        '''
        for r_frac_thresh , r_core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
            r_OTUs , r_biom = r_core_OTUs_biom
            
            ndb_custom_key_r_frac_thres = ndb_custom_key + '~~~~' + r_frac_thresh
            
            if ndb_custom_key_r_frac_thres in local_dict_frac_thresh_otus:
                print "Why do you have same fraction thresholds repeating?"
                #local_dict_frac_thresh_otus[ndb_custom_key_r_frac_thres].append(r_OTUs)
            else:
                local_dict_frac_thresh_otus[ndb_custom_key_r_frac_thres] = r_OTUs

    return local_dict_frac_thresh_otus ## change this to something useful



def shuffle_dicts(a): ## takes dictionary
    keys = a.keys() #If items(), keys(), values() are called with no intervening modifications to the dictionary, the lists will directly correspond.
    values, lengths = get_values_from_dict(a)
    random.shuffle(values)
    new_values = divide_list(values, lengths)
    categ_samples_dict_shuffled = dict(zip(keys, new_values))
    return categ_samples_dict_shuffled


def check_map_file_has_two_groups(a): ## takes iternumb which is number of random iteration and dictionary
    values, lengths = get_values_from_dict(a)
    if len(lengths) > 2:
        return "No"
    return "Yes"


def get_values_from_dict(a):
    values = list()
    lengths = list()
    for i in a.values():
        v = i.split(',')
        values.extend(v)
        lengths.append(len(v)) # sizes of original lists
    return values, lengths 


def compile_results(otus, DELIM): # get unique elements from last column (otus)
    taxon_list = list() # this may be a json or list
    result_ = otus #json.loads(otus)
    for l in result_:
        l = l.strip()
        contents = l.split(DELIM)
        if '#' in l or not l:
            continue
        taxon_list.append(contents[1])
    #print len(taxon_list), len(list(set(taxon_list)))
    return list(set(taxon_list))


def calc_freq_elem_list(a):
    counter=collections.Counter(a)
    return counter


def validate_inputs(ndb_custom_key, user_args, otu_table_biom, c, group, mapping_info_list, p_val_adj, DELIM, NTIMES, OUTPFILE, to_email):
    # find index of SampleID and category to be summarized  # e.g. swg or non-swg
    labels = mapping_info_list[0].split(DELIM)
    indx_sampleid = indx_categ = ''

    errors_list = list()

    try: # this is from the mapping file
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        errors_list.append("'#SampleID' not in the headers of the sample <-> group info file")

    try:
        indx_categ = labels.index(c)
    except ValueError:
        errors_list.append("'%s' not in the headers of the sample <-> group info file" %c)
    
    
    if int(NTIMES) % 50 != 0:
        errors_list.append("Number of randomizations requested is not multiple of 50. Kindly rerun job")


    return (indx_sampleid , indx_categ , errors_list)


# run core microbiome on original data
def run_true_data(OUTPFILE, otu_table_biom, mapping_info_list, c, group, DELIM, ndb_custom_key):
    o_dir = 'true_result' + OUTPFILE
    result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list, c, group) 
    # result dict has following keys 'fractions_for_core' , 'otu_counts' , 'frac_thresh_core_OTUs_biom'

    true_result_frac_thresh_otus_dict = dict() # compile original results. The key is frac_threshold, value is a list of unique otus
    for frac_thresh , core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
        OTUs , biom = core_OTUs_biom # this is a tuple of otus and biom
        true_result_frac_thresh_otus_dict[frac_thresh] = compile_results(OTUs, DELIM)

    Result_TrueDict(parent=ndb.Key(Result_TrueDict, 'fatherresultstrue'), \
                     idx= ndb_custom_key, true_results = true_result_frac_thresh_otus_dict).put()
    print "Processed %s fraction thresholds for true data" % str(len(true_result_frac_thresh_otus_dict))

    return true_result_frac_thresh_otus_dict


def create_shuffled_dicts_no_datastore(i, categ_samples_dict, ndb_custom_key, otu_table_biom):
        '''
        The following section creates shuffled dictionaries and puts them in datastore
        '''

        shuffled_dict = shuffle_dicts(categ_samples_dict)
        ndb_custom_key_entry = ndb_custom_key + '~~~~' + str(i)
        #local_dict_frac_thresh_otus = shuffle_dict_coremic_serial_dict_no_datastore(shuffled_dict, ndb_custom_key, ndb_custom_key_entry, otu_table_biom)
        local_dict_frac_thresh_otus = shuffle_dict_coremic_serial_dict_datastore(shuffled_dict, ndb_custom_key, ndb_custom_key_entry, otu_table_biom)
        
        #RandomDict(parent=ndb.Key(RandomDict, 'father'), idx= ndb_custom_key, entry_id = ndb_custom_key_entry, dict= shuffled_dict).put()
        return local_dict_frac_thresh_otus


def shuffle_dict_coremic_serial_dict_datastore(rand_mapping_info_dict, ndb_custom_key, ndb_custom_key_entry, otu_table_biom):
        '''
        get the randomized dict and run core microbiome
        '''    

        local_dict_frac_thresh_otus = dict()

        OUTPFILE, c, group, rand_iter_numb = ndb_custom_key_entry.split('~~~~') #Zen-outputMon-07-Mar-2016-01:36:13-AM~~~~Plant~~~~Sw~~~~2
        rand_mapping_info_list = convert_shuffled_dict_to_str(rand_mapping_info_dict, c)
        rand_o_dir = rand_iter_numb + OUTPFILE
        result = exec_core_microb_cmd(otu_table_biom, rand_o_dir, rand_mapping_info_list, c, group)

        '''
         arrange the results to look pretty
        '''
        for r_frac_thresh , r_core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
            r_OTUs , r_biom = r_core_OTUs_biom
            
            ndb_custom_key_r_frac_thres = ndb_custom_key + '~~~~' + r_frac_thresh
            
            if ndb_custom_key_r_frac_thres in local_dict_frac_thresh_otus:
                print "Why do you have same fraction thresholds repeating?"
                #local_dict_frac_thresh_otus[ndb_custom_key_r_frac_thres].append(r_OTUs)
            else:
                local_dict_frac_thresh_otus[ndb_custom_key_r_frac_thres] = r_OTUs

        Result_RandomDict(parent=ndb.Key(Result_RandomDict, 'fatherresults'), \
                     idx= ndb_custom_key, otus = local_dict_frac_thresh_otus).put()


        return local_dict_frac_thresh_otus ## change this to something useful


def calc_significance(indx_sampleid , indx_categ , errors_list, otu_table_biom, c, group, mapping_info_list, p_val_adj, DELIM, NTIMES, OUTPFILE, to_email):

    #global ndb_custom_key
    ndb_custom_key = OUTPFILE + '~~~~' + c + '~~~~' + group  # this is to query all entries in this run
    user_args = 'You selected the following parameters:\nFactor: %s\nGroup: %s\nPval correction: %s\n# of randomizations: %s\n\n\n'  %(c, group, p_val_adj, NTIMES)


    categ_samples_dict = list_to_dict(mapping_info_list, DELIM, ',', "current", indx_categ, indx_sampleid)
    
    if check_map_file_has_two_groups(categ_samples_dict) == "No":
        errors_list.append('\nERROR: Following code divides samples in >TWO groups. Change the mapping file to only have two groups (e.g. A vs D)\n')


    if len(errors_list) > 0: # email the error and quit, no point to continue further
        send_results_as_email(ndb_custom_key, user_args, '\n'.join(errors_list), to_email)
        # put code here so that the code doesn't run further

    
    out_str = ''    

    glob_qry_entries_in_result_rand_dict =  dict()
    counter = 1

    print "Creating shuffled dicts"
    for i in range(NTIMES):
        create_shuffled_dicts_no_datastore(i, categ_samples_dict, ndb_custom_key, otu_table_biom)


def compile_all_results_perform_sign_calc(ndb_custom_key, glob_qry_entries_in_result_rand_dict, user_args, to_email, p_val_adj, DELIM, true_result_frac_thresh_otus_dict, NTIMES):
     
    '''
    the following section compiles results from the Result Datatstore and calculates stats.
    '''

    print "Compiling results"
    sign_results = 'Significant results:\nOTU\tFreq. in randomized data\tpval=freq/times randomized\t%s corrected pval\n' %p_val_adj
    p_val = 0.05 


    # compile results; print the number of random occurances for each true core microbiome otu (checks significance)
    for frac_s in [75, 80, 85, 90, 95, 100]:
        sign_results += '\n#Frac thresh %s\n' % str(frac_s)
        
        ndb_custom_key_qury_id = ndb_custom_key + '~~~~' + str(frac_s)
        
        # this number should be equal to the number of randomizations
        # qry_entries_in_result_rand_dict is a list of list, the internal list is tab-delimited OTU# and OTU name
        qry_entries_in_result_rand_dict = glob_qry_entries_in_result_rand_dict[ndb_custom_key_qury_id]
        print 'Number of results from randomized dict for ' , frac_s , '% threshold = ' , len(qry_entries_in_result_rand_dict)
        #print qry_entries_in_result_rand_dict

        # compile the results from randomization
        # this returns a list of list i.e. collects the unique set of core taxa (OTU name) from each randomized data
        taxons_ = list()
        for q in qry_entries_in_result_rand_dict:
            taxons_.append(compile_results(q, DELIM))
        #print taxons_
        
        all_results_otu_list = [item for sublist in taxons_ for item in sublist] #  taxons_ -> sublist -> item
        #print all_results_otu_list

        # calculate freq of otu being a core microb from the randomizations
        randomized_otus = calc_freq_elem_list(all_results_otu_list) # dict of otus and freq of occurance from random data

        # check significance
        signif_core_microb_otu_dict = collections.OrderedDict()
        for o in true_result_frac_thresh_otus_dict[str(frac_s)]:
            #print o
            freq = 0.1
            if o in randomized_otus: # the else for this means it was not observed even once in the randomized data!
                freq = int(randomized_otus[o])
            otus_pval = freq/float(NTIMES)
            if otus_pval < p_val:
                otu = '%s\t%s\t%s\n' % (o, freq, otus_pval)
                #print otu
                #sign_results += otu
                signif_core_microb_otu_dict[o] = otus_pval

        # check if there is at least one significant entry so far:
        if len(signif_core_microb_otu_dict) == 0:
            continue # go to next iteration
        else:
            print "There are %s signif core microb without corrections"  % str(len(signif_core_microb_otu_dict))

        #print signif_core_microb_otu_dict.values()
        # adjust the pvalues of significant otus for multiple testing
        new_p_vals = list()
        if p_val_adj == 'bf':
            new_p_vals = correct_pvalues_for_multiple_testing(signif_core_microb_otu_dict.values(), "Bonferroni")
        elif p_val_adj == 'bh':
            new_p_vals = correct_pvalues_for_multiple_testing(signif_core_microb_otu_dict.values(), "Benjamini-Hochberg")
        elif p_val_adj == 'none':
            new_p_vals = signif_core_microb_otu_dict.values()

        
        counter = 0
        for o in signif_core_microb_otu_dict.keys():
            freq = int(randomized_otus[o])
            otus_pval = freq/float(NTIMES) # p value before correction (from randomized runs)
            new_p_v = new_p_vals[counter]  # p value after correction 
            if new_p_v < p_val:
                otu = '%s\t%s\t%s\t%s\n' % (o, freq, otus_pval, new_p_v)
                sign_results += otu
            counter += 1

    
    #print sign_results
    send_results_as_email(ndb_custom_key, user_args, sign_results, to_email)
    return sign_results #r_out_str #result['otu_counts'] # '1'


# the user id needs to be changed to that input by the user
def send_results_as_email(timestmp, user_args, msg, to_email):
    subj = "Your data from %s has been processed" %timestmp
    message = mail.EmailMessage(sender="Core microbiome Support <richieangel@gmail.com>",
                            subject=subj)

    email_to_ =  "User <" + to_email + ">" #"User <richrr@vt.edu>" #Albert Johnson <Albert.Johnson@example.com>
    message.to = email_to_

    msg_str = """
Dear User:

Your data has been processed. Thanks for using this tool.

Please email us if you have any questions.

The Core Microbiome Team

"""
    
    msg_str += user_args
    
    msg_str += msg

    message.body = msg_str

    message.send()


# http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
# the pvalues do not have to be sorted, their order is maintained in the results
def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty                                                                        
    pvalues = array(pvalues) 
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":                                                                   
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    elif correction_type == "None":
        new_pvalues = 1 * pvalues                                                                                         
    return new_pvalues


def convert_shuffled_dict_to_str(DICT, categ):
    file_str_list = list()
    file_str = "%s\t%s\n" % ('#SampleID' , categ)
    file_str_list.append(file_str)

    for k, v in DICT.items():
        for i in v:
            f_str = "%s\t%s\n" % (i , k)
            file_str_list.append(f_str)

    return file_str_list


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))



MAIN_PAGE_HTML = """\
<html>
  <body>
    <form action="{0}" enctype="multipart/form-data" method="post">
    <b> 'Calculate significance of core microbiome' </b>
	<p>
		Type some text to use for your output (if you like):<br>
		<input type="text" name="content" value="Zen-output" size="30">
	</p>
	<p>
		Please specify the BIOM file containing an OTU table:<br>
		make sure the header "taxonomy" is present<br>
		<input type="file" name="datafile" size="40">
	</p>
	<p>
		Please specify the <b>tab-delimited</b> file containing sample <-> group info:<br>
		e.g. the group of interest "Person" has to be binary ("Good" vs "Bad")<br>
		<input type="file" name="groupfile" size="40">
	</p>
	<p>
		Enter EXACT string of the interest group on which you want core microbiome analysis to be performed:<br>
		e.g., To calculate core microbiome of "Good" under the "Person" group, enter "Person:Good"<br>
		<input type="text" name="group" value="Person:Good" size="30">
	</p>
	<p>
		Enter number of randomizations to be performed (USE multiples of 50; 500 to 1,500):<br>
		<input type="text" name="random" value="1000" size="30">
	</p>
	<p>
	  <input type="radio" name="pvaladjmethod" value="bf"> Bonferroni<br>
	  <input type="radio" name="pvaladjmethod" value="bh" checked> Benjamini Hochberg<br>
	  <input type="radio" name="pvaladjmethod" value="none"> None
	</p>
	<p>
		Enter the email address where you want your results emailed:<br>
		Your results are NOT stored, you will only get these results via email<br>
		<input type="text" name="email" value="user@domain.com" size="50">
	</p>
	<div>
		<input type="submit" value="Send">
	</div>
    </form>
    
    <div>
		Developed by Rodrigues RR and Williams MA, Virginia Tech. <br>
		Acknowledgements: QIIME and GAE. <br>
		This tool is free to use for non-profit or research purposes. <br>
    See here for help section. <br>
	</div>

  </body>
</html>
"""


class MainPage(webapp2.RequestHandler):
    def get(self):
        upload_url = '/sign'
        self.response.write(MAIN_PAGE_HTML.format(upload_url))


#class Guestbook(blobstore_handlers.BlobstoreUploadHandler,webapp2.RequestHandler):
class Guestbook(webapp2.RequestHandler):
    def post(self):
        #otu_table_biom = self.get_uploads()[0]
        #group_info = self.get_uploads()[1]
        otu_table_biom = self.request.get('datafile') #automatically reads file
    	group_info = self.request.get('groupfile') #automatically reads file
        
        #print self.get_uploads()
        
        #global DELIM, NTIMES, OUTPFILE
        DELIM = '\t'
    	factor, group = self.request.get('group').split(':')
    	NTIMES = int(self.request.get('random'))
    	datetime = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())
    	
    	## added a random alpha numeric to avoid conflict with another (simultaneous) user request.
        OUTPFILE = self.request.get('content') + datetime + id_generator()

    	p_val_adj = self.request.get('pvaladjmethod')

        to_email = self.request.get('email')

        ndb_custom_key_o = OUTPFILE + '~~~~' + factor + '~~~~' + group  # this is to query all entries in this run

        # make a dict and insert in ndb as a json property
        local_string_hack_dict = { "otu_table_biom_key" : ndb_custom_key_o, 
                     "g_info_not_list" : group_info, 
                     "factor" : factor, 
                     "group" : group, 
                     "p_val_adj" : p_val_adj, 
                     "delim" : DELIM, 
                     "ntimes" : str(NTIMES), 
                     "outpfile" : OUTPFILE, 
                     "to_email" : to_email }
        
        ''' 
        temporary hack to clear out the Datastore     http://stackoverflow.com/questions/1062540/how-to-delete-all-datastore-in-google-app-engine
        '''
        #ndb.delete_multi(RandomDict.query().fetch(keys_only=True)) 
        ndb.delete_multi(Result_RandomDict.query().fetch(keys_only=True)) 
        #ndb.delete_multi(ResultFile.query().fetch(keys_only=True)) 
        ndb.delete_multi(OriginalBiom.query().fetch(keys_only=True)) 
        ndb.delete_multi(Result_TrueDict.query().fetch(keys_only=True)) 


        OriginalBiom(id='origbiom').put()  # the datastore of original biom
        OriginalBiom(parent=ndb.Key(OriginalBiom, 'origbiom'), 
            idx= ndb_custom_key_o, biom = otu_table_biom, params_str=local_string_hack_dict).put()

        Result_RandomDict(id='fatherresults').put() # the datastore of results from random dicts

        Result_TrueDict(id='fatherresultstrue').put() # the datastore of results from random dicts

        numb_tasks = int(NTIMES)/50
        # try to break the randomizations into n tasks of 50 randomizations each and then run them one by one
        # finally run the significance calculation
        # http://stackoverflow.com/questions/4224564/calling-a-script-after-tasks-queue-is-empty 
        taskqueue.add(url="/process_data", params={'otu_table_biom_key': ndb_custom_key_o , 'true_random' : 'true'},
                 retry_options=TaskRetryOptions(task_retry_limit=0, task_age_limit=1),
                countdown=1)
        for i in range(numb_tasks):
            taskqueue.add(url="/process_data", params={'otu_table_biom_key': ndb_custom_key_o , 'true_random' : 'random'},
                 retry_options=TaskRetryOptions(task_retry_limit=0, task_age_limit=1),
                countdown=1)

        taskqueue.add(url="/process_results", params={'otu_table_biom_key': ndb_custom_key_o, 'numb_tasks' : numb_tasks},
                retry_options=TaskRetryOptions(task_retry_limit=0, task_age_limit=1),
                countdown=1)
        
        self.redirect('/')
        

class ProcessData(webapp2.RequestHandler):
    def post(self):
        otu_table_biom_o = self.request.get("otu_table_biom_key")
        true_random = self.request.get("true_random")
        '''
    def get(self, photo_key):
        '''
        
        ## get rid of this orig biom and use only the method
        ## also move the true results to a separate task to avoid redoing work
        user_args, to_email, p_val_adj, DELIM, NTIMES, otu_table_biom, g_info_list, factor, group, OUTPFILE = get_required_params_from_orig_dict(otu_table_biom_o)
            
        indx_sampleid , indx_categ , errors_list = validate_inputs("ndb_custom_key", "user_args", otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)
        if len(errors_list) > 0: # temp hack since blobstore randomly swaps file order during upload
            tmp = g_info_list
            g_info_list = otu_table_biom
            otu_table_biom = tmp
            # retry with swapped files
            indx_sampleid , indx_categ , errors_list = validate_inputs("ndb_custom_key", "user_args", otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)
            if len(errors_list) > 0: # just give up on this
                print '\n'.join(errors_list)
                # put code here so that the code doesn't run further
                sys.exit(0)
            else:
                print 'Swapping files worked!'
        else:
            print 'No file swapping needed!'

        #print 'OTU' , otu_table_biom, '\n', 'GINFO', g_info_list
        ndb_custom_key = OUTPFILE + '~~~~' + factor + '~~~~' + group  # this is to query all entries in this run

        if true_random == "true":
            run_true_data(OUTPFILE, otu_table_biom, g_info_list, factor, group, DELIM, ndb_custom_key)
        elif true_random == "random":
            calc_significance(indx_sampleid , indx_categ , errors_list, otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)


def get_required_params_from_orig_dict(otu_table_biom_o):
        qry_entries_in_origbiom = OriginalBiom.query(OriginalBiom.idx == otu_table_biom_o, ancestor=ndb.Key(OriginalBiom, 'origbiom'))
        
        if int(qry_entries_in_origbiom.count()) == 1:
            print "Read single entry from OriginalBiom datastore!"
        else:
            # do something useful here
            print "Single entry not found from OriginalBiom datastore, some error!"
        
        for q in qry_entries_in_origbiom:
            print "Running 1st attempt"
            q_dict = q.to_dict()
           
            params = q_dict['params_str']  # this is a dictionary
            factor = params["factor"]
            group = params["group"]
            p_val_adj = params["p_val_adj"]
            DELIM = params["delim"]
            NTIMES = str(50) #params["ntimes"])
            OUTPFILE = params["outpfile"]
            to_email = params["to_email"]
            
            otu_table_biom = q_dict['biom']        
             
            g_info_list = params["g_info_not_list"].split('\n')
            
            user_args = 'You selected the following parameters:\nFactor: %s\nGroup: %s\nPval correction: %s'  %(factor, group, p_val_adj)

            return user_args, to_email, p_val_adj, DELIM, NTIMES, otu_table_biom, g_info_list, factor, group, OUTPFILE


class ProcessResults(webapp2.RequestHandler):
    def post(self):
        otu_table_biom_o = self.request.get("otu_table_biom_key")
        numb_tasks = self.request.get("numb_tasks")

        true_result_frac_thresh_otus_dict = dict()
        qry_entries_in_result_Truedict = Result_TrueDict.query(Result_TrueDict.idx == otu_table_biom_o, ancestor=ndb.Key(Result_TrueDict, 'fatherresultstrue'))
        
        if int(qry_entries_in_result_Truedict.count()) == 1:
            print "Read single entry from Truedict datastore!"
        else:
            # do something useful here
            print "Single entry not found from Truedict datastore, some error!"
        
        for t in qry_entries_in_result_Truedict:
            t_dict = t.to_dict()
            true_result_frac_thresh_otus_dict = t_dict['true_results']
            break

        qry_entries_in_result_randomdict = Result_RandomDict.query(Result_RandomDict.idx == otu_table_biom_o, ancestor=ndb.Key(Result_RandomDict, 'fatherresults'))
        qry_entries_in_result_randomdict_count = Result_RandomDict.query(Result_RandomDict.idx == otu_table_biom_o, ancestor=ndb.Key(Result_RandomDict, 'fatherresults')).count(limit=100000)
        print "Counts" , qry_entries_in_result_randomdict_count, qry_entries_in_result_Truedict.count(), qry_entries_in_result_randomdict.count()
        
        #print true_result_frac_thresh_otus_dict
        
        if int(qry_entries_in_result_randomdict_count) == (int(numb_tasks)*50):
            print "Previous work completed, can move for final stage!"
            # merge all the available dictionaries into one
            glob_qry_entries_in_result_rand_dict =  dict()

            for q in qry_entries_in_result_randomdict:
                q_dict = q.to_dict()
                local_dict_frac_thresh_otus = q_dict['otus']
                for ndb_custom_key_r_frac_thres , r_OTUs in local_dict_frac_thresh_otus.items():
                    if ndb_custom_key_r_frac_thres in glob_qry_entries_in_result_rand_dict:
                        glob_qry_entries_in_result_rand_dict[ndb_custom_key_r_frac_thres].append(r_OTUs)
                    else:
                        glob_qry_entries_in_result_rand_dict[ndb_custom_key_r_frac_thres] = [r_OTUs]
            tmp_user_args, to_email, p_val_adj, DELIM, NTIMES, otu_table_biom, g_info_list, factor, group, OUTPFILE = get_required_params_from_orig_dict(otu_table_biom_o)
            NTIMES = int(numb_tasks)*50
            
            user_args = tmp_user_args + "\n# of randomizations: " + str(NTIMES) + "\n\n\n"
            compile_all_results_perform_sign_calc(otu_table_biom_o, glob_qry_entries_in_result_rand_dict, user_args, to_email, p_val_adj, DELIM, true_result_frac_thresh_otus_dict, NTIMES)  
            # may want to purge remaining tasks, be careful since you do not want to delete someone
            # else's tasks
            # Purge entire queue...
            purgeq = taskqueue.Queue('default')
            purgeq.purge()


        else:
            # do something useful here
            print "Waiting for previous tasks to finish!"
            #time.sleep(55) # allow some time for the other process to finish
            taskqueue.add(url="/process_results", params={'otu_table_biom_key': otu_table_biom_o, 'numb_tasks' : numb_tasks},
                 retry_options=TaskRetryOptions(task_retry_limit=0, task_age_limit=1),
                countdown=1)

      



app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/sign', Guestbook),
    ('/process_data', ProcessData),
    ('/process_results', ProcessResults),
], debug=True)


