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

from google.appengine.ext.blobstore import BlobKey

#from google.appengine.ext import deferred

# to do the parallelism: async, futures or mapreduce
# https://cloud.google.com/appengine/docs/python/datastore/async#Working_with_the_Async_Datastore_API
# futures is for python 3, so cannot be used for my 2.7
# mapreduce should work; if it doesn't, try the async datastore api


#http://stackoverflow.com/questions/11849456/how-to-filter-datastore-data-before-mapping-to-cloud-storage-using-the-mapreduce
#http://stackoverflow.com/questions/23508116/appengine-mapreduce-how-to-filter-structuredproperty-while-using-datastore-input

# This datastore model keeps track of which users uploaded which photos.
class UserPhoto(ndb.Model):
    request_id = ndb.StringProperty()
    otu_table_biom_blob_key = ndb.BlobKeyProperty()
    group_info_blob_key = ndb.BlobKeyProperty()
    params_str = ndb.JsonProperty()


#"filters": [("idx", "=", ndb_custom_key)]
class ShuffleDictPipeline(base_handler.PipelineBase):
  def run(self, ndb_custom_key, otu_table_biom):
    output = yield mapreduce_pipeline.MapperPipeline(
      "calc_shuff_core_microb",
      "coremic.shuffle_dict_coremic_map",
      "mapreduce.input_readers.DatastoreInputReader", 
      output_writer_spec="mapreduce.output_writers.GoogleCloudStorageConsistentOutputWriter",
      params={
        "input_reader":{
          "entity_kind":  "coremic.RandomDict",
          "ndb_custom_key" : ndb_custom_key,
          "otu_table_biom" : otu_table_biom
        },
        "output_writer":{
          "filesystem": "gs",
          "bucket_name": "coremicrobucket"
        }
      },
      shards=250) #200
    #print output
    #with pipeline.After(output):
    yield CloudStorageWriter(output)


class ResultFile(ndb.Model):
  file_name = ndb.StringProperty()
  date = ndb.DateTimeProperty(auto_now_add=True)


## this is useless information for now, fix it or delete it
class CloudStorageWriter(base_handler.PipelineBase):
    def run(self, csv_output):
      # Store all the file names
      files = [str(f.replace('/gs/', 'gs://')) for f in csv_output]
      for f in files:
        entry = ResultFile(file_name=f)
        entry.put()


# the actual key is automatically generated
class RandomDict(ndb.Model):
  idx = ndb.StringProperty() # the run id
  entry_id = ndb.StringProperty() # the entry along with the run id
  dict = ndb.JsonProperty() # the (shuffled) dictionary as a json


## for every random dict entry, it has different thresholds
# the actual key is automatically generated
class Result_RandomDict(ndb.Model):
  idx = ndb.StringProperty() # the run id
  frac_thresh = ndb.StringProperty() # the frac threshold along with the run id
  entry_id = ndb.StringProperty() # the frac threshold and entry along with the run id
  otus = ndb.JsonProperty() # the core otus from the (shuffled) dictionary as a json
  biom = ndb.JsonProperty() # the core biom from the (shuffled) dictionary as a json


# the actual key is automatically generated
class OriginalBiom(ndb.Model):
  idx = ndb.StringProperty() # the run id
  biom = ndb.JsonProperty() # the core biom from the original input as a json


# this cane be changed later to allow splitting in more than two groups
def divide_list(a, lengths):
    return a[:lengths[0]] , a[lengths[0]:]  # a[start:end] # items start through end-1


def shuffle_dict_coremic_map(entity):
    '''
     get the params and any information that was passed
    '''
    ctx = context.get()
    params = ctx.mapreduce_spec.mapper.params
    #print params #{u'input_reader': {u'ndb_custom_key': u'Zen-outputTue-05-Apr-2016-11:31:08-PM43MDJ2~~~~Plant~~~~Sw', u'entity_kind': u'coremic.RandomDict', u'otu_table_biom' : u'too big dict...to show here'}, u'output_writer': {u'filesystem': u'gs', u'bucket_name': u'coremicrobucket'}}
    ndb_custom_key = params['input_reader']["ndb_custom_key"]
    otu_table_biom = params['input_reader']["otu_table_biom"]
    #print otu_table_biom

    '''
     get the randomized dict and run core microbiome
    '''    
    entty = entity.to_dict()
    #print ndb_custom_key , 'aaaaaaaaaaaaaaa'
    #print entty['idx'] , 'AAAAAAAAAAAAAAAA'
    r_out_str = ''
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
            #r_out_str += "Frac Threshold %s:\n%s\n%s\n\n" % (r_frac_thresh, ''.join(r_OTUs), r_biom)
            #r_out_str += "Frac Threshold %s:\n%s\n\n" % (r_frac_thresh, ''.join(r_OTUs))
            
            ndb_custom_key_r_frac_thres = ndb_custom_key + '~~~~' + r_frac_thresh
            ndb_custom_key_r_frac_thres_entry = ndb_custom_key_r_frac_thres + '~~~~' + rand_iter_numb

            Result_RandomDict(parent=ndb.Key(Result_RandomDict, 'fatherresults'), \
                     idx= ndb_custom_key, frac_thresh = ndb_custom_key_r_frac_thres, \
                     entry_id = ndb_custom_key_r_frac_thres_entry, \
                     otus = r_OTUs, biom = r_biom).put()

    yield '1' ## change this to something useful


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

    return (indx_sampleid , indx_categ , errors_list)


# run core microbiome on original data
def run_true_data(OUTPFILE, otu_table_biom, mapping_info_list, c, group, DELIM):
    o_dir = 'true_result' + OUTPFILE
    result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list, c, group) 
    # result dict has following keys 'fractions_for_core' , 'otu_counts' , 'frac_thresh_core_OTUs_biom'

    true_result_frac_thresh_otus_dict = dict() # compile original results. The key is frac_threshold, value is a list of unique otus
    for frac_thresh , core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
        OTUs , biom = core_OTUs_biom # this is a tuple of otus and biom
        true_result_frac_thresh_otus_dict[frac_thresh] = compile_results(OTUs, DELIM)
    return true_result_frac_thresh_otus_dict


def create_shuffled_dicts(NTIMES, categ_samples_dict, ndb_custom_key):
    '''
    The following section creates shuffled dictionaries and puts them in datastore
    '''
    # do this to get consistency in queries, else the number of returned results vary
    #http://stackoverflow.com/questions/14630886/ndb-and-consistency-why-is-happening-this-behavior-in-a-query-without-a-parent
    RandomDict(id='father').put()  # the datastore of random dicts

    for i in range(NTIMES):
        shuffled_dict = shuffle_dicts(categ_samples_dict)
        ndb_custom_key_entry = ndb_custom_key + '~~~~' + str(i)
        RandomDict(parent=ndb.Key(RandomDict, 'father'), idx= ndb_custom_key, entry_id = ndb_custom_key_entry, dict= shuffled_dict).put()
    return 1


def calc_significance(indx_sampleid , indx_categ , errors_list, otu_table_biom, c, group, mapping_info_list, p_val_adj, DELIM, NTIMES, OUTPFILE, to_email):

    #global ndb_custom_key
    ndb_custom_key = OUTPFILE + '~~~~' + c + '~~~~' + group  # this is to query all entries in this run
    user_args = 'You selected the following parameters:\nFactor: %s\nGroup: %s\nPval correction: %s\n# of randomizations: %s\n\n\n'  %(c, group, p_val_adj, NTIMES)


    categ_samples_dict = list_to_dict(mapping_info_list, DELIM, ',', "current", indx_categ, indx_sampleid)
    
    if check_map_file_has_two_groups(categ_samples_dict) == "No":
        errors_list.append('\nERROR: Following code divides samples in >TWO groups. Change the mapping file to only have two groups (e.g. A vs D)\n')


    if len(errors_list) > 0: # email the error and quit, no point to continue further
        send_results_as_email(ndb_custom_key, user_args, '\n'.join(errors_list), to_email)

    
    out_str = ''
    true_result_frac_thresh_otus_dict = run_true_data(OUTPFILE, otu_table_biom, mapping_info_list, c, group, DELIM)
    
    ''' 
    temporary hack to clear out the Datastore     http://stackoverflow.com/questions/1062540/how-to-delete-all-datastore-in-google-app-engine
    '''
    ndb.delete_multi(RandomDict.query().fetch(keys_only=True)) 
    ndb.delete_multi(Result_RandomDict.query().fetch(keys_only=True)) 
    ndb.delete_multi(ResultFile.query().fetch(keys_only=True)) 
    ndb.delete_multi(OriginalBiom.query().fetch(keys_only=True)) 


    create_shuffled_dicts(NTIMES, categ_samples_dict, ndb_custom_key)


    '''
    The following section calculates core microbiome for each shuffled dictionaries and puts the
    results in the result datastore
    '''
 
    Result_RandomDict(id='fatherresults').put() # the datastore of results from random dicts
    # query entries with same ndb_custom_key
    #qry_entries_in_rand_dict = RandomDict.query(RandomDict.idx == ndb_custom_key, ancestor=ndb.Key(RandomDict, 'father'))  
    #print qry_entries_in_rand_dict.count()

    ## using mapreduce to parallelize the core microbiome on random dicts
    #shuffled_core_mic = ShuffleDictPipeline()#ndb_custom_key, otu_table_biom) 
    shuffled_core_mic = ShuffleDictPipeline(ndb_custom_key, otu_table_biom) 
    shuffled_core_mic.start()


    # this keeps the below code from running until mapreduce is finished    
    #time.sleep(55)
    
    
    # Later on, see if it's done.
    my_pipeline = shuffled_core_mic.pipeline_id
    shuffled_core_mic = ShuffleDictPipeline.from_id(my_pipeline)
    status_var = shuffled_core_mic.has_finalized
    #print status_var
    while not status_var:
        shuffled_core_mic = ShuffleDictPipeline.from_id(my_pipeline)
        status_var = shuffled_core_mic.has_finalized
        print status_var
        if status_var:
            print "........Breaking........"
            break
        print "........Waiting.........sleep 30 secs"
        time.sleep(30)

    print "........Done Waiting........."

     
    '''
    the following section compiles results from the Result Datatstore and calculates stats.
    '''

    sign_results = 'Significant results:\nOTU\tFreq. in randomized data\tpval=freq/times randomized\t%s corrected pval\n' %p_val_adj
    p_val = 0.05 

    # compile results; print the number of random occurances for each true core microbiome otu (checks significance)
    for frac_s in [75, 80, 85, 90, 95, 100]:#for frac_s in [100]:
        sign_results += '#Frac thresh %s\n' % str(frac_s)
        
        ndb_custom_key_qury_id = ndb_custom_key + '~~~~' + str(frac_s)
       
        qry_entries_in_result_rand_dict = Result_RandomDict.query(Result_RandomDict.frac_thresh == ndb_custom_key_qury_id, ancestor=ndb.Key(Result_RandomDict, 'fatherresults'))  
        #print 'xxxxxxxxxxxxx' , qry_entries_in_result_rand_dict.count()

        # compile the results from randomization
        # this returns a list of list i.e. collects the unique set of core taxa from each randomized data
        taxons_ = list()
        for q in qry_entries_in_result_rand_dict:
            res = q.to_dict()
            taxons_.append(compile_results(res['otus'], DELIM))
            
        all_results_otu_list = [item for sublist in taxons_ for item in sublist] #  taxons_ -> sublist -> item

        # calculate freq of otu being a core microb from the randomizations
        randomized_otus = calc_freq_elem_list(all_results_otu_list) # dict of otus and freq of occurance from random data

        # check significance
        signif_core_microb_otu_dict = collections.OrderedDict()
        for o in true_result_frac_thresh_otus_dict[str(frac_s)]:
            if o in randomized_otus:
                freq = int(randomized_otus[o])
                otus_pval = freq/float(NTIMES)
                if otus_pval < p_val:
                    #otu = '%s\t%s\t%s\n' % (o, freq, otus_pval)
                    #sign_results += otu
                    signif_core_microb_otu_dict[o] = otus_pval

        # check if there is at least one significant entry so far:
        if len(signif_core_microb_otu_dict) == 0:
            continue # go to next iteration
        else:
            pass    

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

    
    print sign_results
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
		Enter number of randomizations to be performed (500 to 10,000):<br>
		<input type="text" name="random" value="1000" size="30">
	</p>
	<p>
	  <input type="radio" name="pvaladjmethod" value="bf" checked> Bonferroni<br>
	  <input type="radio" name="pvaladjmethod" value="bh"> Benjamini Hochberg<br>
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
        upload_url = blobstore.create_upload_url('/sign')
        self.response.write(MAIN_PAGE_HTML.format(upload_url))


class Guestbook(blobstore_handlers.BlobstoreUploadHandler,webapp2.RequestHandler):
    def post(self):
        otu_table_biom = self.get_uploads()[0]
        group_info = self.get_uploads()[1]
        
        print self.get_uploads()
        
        #global DELIM, NTIMES, OUTPFILE
        DELIM = '\t'
    	factor, group = self.request.get('group').split(':')
    	NTIMES = int(self.request.get('random'))
    	datetime = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())
    	
    	## added a random alpha numeric to avoid conflict with another (simultaneous) user request.
        OUTPFILE = self.request.get('content') + datetime + id_generator()

    	p_val_adj = self.request.get('pvaladjmethod')

        to_email = self.request.get('email')

        # make a dict and insert in ndb as a json property
        local_string_hack_dict = { "otu_table_biom_key" : str(otu_table_biom.key()), 
                     "group_key" : str(group_info.key()), 
                     "factor" : factor, 
                     "group" : group, 
                     "p_val_adj" : p_val_adj, 
                     "delim" : DELIM, 
                     "ntimes" : str(NTIMES), 
                     "outpfile" : OUTPFILE, 
                     "to_email" : to_email }
            
        UserPhoto(id='UserPhoto').put()  # the datastore of original biom
        UserPhoto(parent=ndb.Key(UserPhoto, 'UserPhoto'), 
                request_id= OUTPFILE,
                otu_table_biom_blob_key=otu_table_biom.key(),
                group_info_blob_key=group_info.key(),
                params_str=local_string_hack_dict).put()

        self.redirect('/process_data/%s' % otu_table_biom.key())

class ProcessData(blobstore_handlers.BlobstoreDownloadHandler):
    def get(self, photo_key):
        qry_entries_in_origbiom = UserPhoto.query(UserPhoto.otu_table_biom_blob_key == BlobKey(photo_key), ancestor=ndb.Key(UserPhoto, 'UserPhoto'))
        
        if int(qry_entries_in_origbiom.count()) == 1:
            pass
        else:
            # do something useful here
            pass
        
        for q in qry_entries_in_origbiom:
            q_dict = q.to_dict()
            
            params = q_dict['params_str']  # this is a dictionary
            factor = params["factor"]
            group = params["group"]
            p_val_adj = params["p_val_adj"]
            DELIM = params["delim"]
            NTIMES = str(params["ntimes"])
            OUTPFILE = params["outpfile"]
            to_email = params["to_email"]
            
            otu_table_biom_key = q_dict['otu_table_biom_blob_key']
            group_info_key = q_dict['group_info_blob_key']
            #print otu_table_biom_key , '\n', group_info_key
            
            # read contents of the blob
            blob_reader = blobstore.BlobReader(otu_table_biom_key)
            otu_table_biom = blob_reader.read()
            blob_reader.close()
        
            # read contents of the blob
            blob_reader = blobstore.BlobReader(group_info_key)
            g_info_list = blob_reader.readlines()
            blob_reader.close()
            
            indx_sampleid , indx_categ , errors_list = validate_inputs("ndb_custom_key", "user_args", otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)
            if len(errors_list) > 0: # temp hack since blobstore randomly swaps file order during upload
                tmp = g_info_list
                g_info_list = otu_table_biom
                otu_table_biom = tmp
                # retry with swapped files
                indx_sampleid , indx_categ , errors_list = validate_inputs("ndb_custom_key", "user_args", otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)
                if len(errors_list) > 0: # just give up on this
                    print '\n'.join(errors_list)
                    sys.exit(0)
                else:
                    print 'Swapping files worked!'
            else:
                print 'No file swapping needed!'

            #print g_info_list

            '''
            for i in [otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, NTIMES, OUTPFILE, to_email]:
                print i
            '''
            calc_significance(indx_sampleid , indx_categ , errors_list, otu_table_biom, factor, group, g_info_list, p_val_adj, DELIM, int(NTIMES), OUTPFILE, to_email)
            #self.send_blob(photo_key)


        



app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/sign', Guestbook),
    ('/process_data/([^/]+)?', ProcessData),
], debug=True)


