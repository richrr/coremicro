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
from google.appengine.api import taskqueue
'''
from mapreduce import base_handler

from mapreduce import mapreduce_pipeline
from mapreduce import operation as op
from mapreduce import shuffler

'''
entity_kind
the datastore kind to map over.
keys_only
use a keys_only query.
batch_size
the number of entities to read from the datastore with each batch get.
key_range
a range of keys to return from your query
filters
Any filters to apply to the datastore query.
'''

# to do the parallelism: async, futures or mapreduce
# https://cloud.google.com/appengine/docs/python/datastore/async#Working_with_the_Async_Datastore_API
# futures is for python 3, so cannot be used for my 2.7
# mapreduce should work; if it doesn't, try the async datastore api


#http://stackoverflow.com/questions/11849456/how-to-filter-datastore-data-before-mapping-to-cloud-storage-using-the-mapreduce
#http://stackoverflow.com/questions/23508116/appengine-mapreduce-how-to-filter-structuredproperty-while-using-datastore-input

class ShuffleDictPipeline(base_handler.PipelineBase):
    """ Count the number of occurrences of a character in a set of strings. """

    def run(self, *args, **kwargs):
        """ run """
        ## this is for the input reader and output writer, not about the methods handling them
        ## change them after you change the input reader and output writer
        mapper_params = {
            "entity_kind": "coremic.RandomDict",
            "batch_size": 500,
            "filters": [("idx", "=", ndb_custom_key)]
        }
        reducer_params = {
            "mime_type": "text/plain"
        }
        output = yield mapreduce_pipeline.MapreducePipeline(
            "calc_shuff_core_microb",
            mapper_spec="coremic.shuffle_dict_coremic_map",
            mapper_params=mapper_params,
            reducer_spec="coremic.shuffle_dict_dict_coremic_reduce", 
            reducer_params=reducer_params,
            input_reader_spec="mapreduce.input_readers.DatastoreInputReader", 
            output_writer_spec="mapreduce.output_writers.BlobstoreOutputWriter",  ### this needs to be changed
            shards=16)

        yield StoreOutput(output)


class StoreOutput(base_handler.PipelineBase):
    """ A pipeline to store the result of the MapReduce job in the database. """

    def run(self, output):
        """ run """
        counter = CharacterCounter(count_link=output[0])
        counter.put()


class CharacterCounter(ndb.Model):
    """ A simple model to sotre the link to the blob storing our MapReduce output. """
    count_link = ndb.StringProperty(required=True)


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


# this cane be changed later to allow splitting in more than two groups
def divide_list(a, lengths):
    return a[:lengths[0]] , a[lengths[0]:]  # a[start:end] # items start through end-1


def shuffle_dict_coremic_map(entity):
    entty = entity.to_dict()
    ### filter the queries as per the run id
    ### split the entry id to get the required information.
    #o_dir =  + OUTPFILE
    #result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list, c, group) 
    yield (entty['dict'],1)

   
def shuffle_dict_coremic_reduce(key, val):
    yield (key, val)


def shuffle_dicts(iternumb, a): ## takes iternumb which is number of random iteration and dictionary
    keys = a.keys() #If items(), keys(), values() are called with no intervening modifications to the dictionary, the lists will directly correspond.
    values = list()
    lengths = list()
    for i in a.values():
        v = i.split(',')
        values.extend(v)
        lengths.append(len(v)) # sizes of original lists
    random.shuffle(values)
    if len(lengths) > 2:
        return '\nFollowing code divides samples in TWO groups. Change the mapping file to only have two groups (e.g. A vs D)\n'
    new_values = divide_list(values, lengths)
    categ_samples_dict_shuffled = dict(zip(keys, new_values))
    return categ_samples_dict_shuffled


def compile_results(otus): # get unique elements from last column (otus)
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


def calc_significance(otu_table_biom, c, group, mapping_info_list, p_val_adj):
    # find index of SampleID and category to be summarized  # e.g. swg or non-swg
    labels = mapping_info_list[0].split(DELIM)
    indx_sampleid = indx_categ = ''

    try: # this is from the mapping file
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        return "'#SampleID' not in the headers of the sample <-> group info file"

    try:
        indx_categ = labels.index(c)
    except ValueError:
        return "'%s' not in the headers of the sample <-> group info file" %c

    # run core microbiome on original data
    o_dir = 'true_result' + OUTPFILE
    result = exec_core_microb_cmd(otu_table_biom, o_dir, mapping_info_list, c, group) 
    # Access the result dict with the following keys 'fractions_for_core' , 'otu_counts' , 'frac_thresh_core_OTUs_biom'

    out_str = ''
    true_result_frac_thresh_otus_dict = dict() # compile original
    for frac_thresh , core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
        OTUs , biom = core_OTUs_biom
        out_str += "Frac Threshold %s:\n%s\n%s\n\n" % (frac_thresh, ''.join(OTUs), biom)
        #print frac_thresh
        true_result_frac_thresh_otus_dict[frac_thresh] = compile_results(OTUs)

    categ_samples_dict = list_to_dict(mapping_info_list, DELIM, ',', "current", indx_categ, indx_sampleid)
   
    ''' 
    temporary hack to clear out the Datastore
    http://stackoverflow.com/questions/1062540/how-to-delete-all-datastore-in-google-app-engine
    '''
    ndb.delete_multi(RandomDict.query().fetch(keys_only=True)) 
    ndb.delete_multi(Result_RandomDict.query().fetch(keys_only=True)) 

    '''
    The following section creates shuffled dictionaries and puts them in datastore
    '''
    # do this to get consistency in queries, else the number of returned results vary
    #http://stackoverflow.com/questions/14630886/ndb-and-consistency-why-is-happening-this-behavior-in-a-query-without-a-parent
    RandomDict(id='father').put()  # the datastore of random dicts

    global ndb_custom_key
    ndb_custom_key = OUTPFILE + '~~~~' + c + '~~~~' + group  # this is to query all entries in this run
    for i in range(NTIMES):
        shuffled_dict = shuffle_dicts(i, categ_samples_dict)
        ndb_custom_key_entry = ndb_custom_key + '~~~~' + str(i)
        RandomDict(parent=ndb.Key(RandomDict, 'father'), idx= ndb_custom_key, entry_id = ndb_custom_key_entry, dict= shuffled_dict).put()


    '''
    The following section calculates core microbiome for each shuffled dictionaries and puts the
    results in the result datastore
    '''
    Result_RandomDict(id='fatherresults').put() # the datastore of results from random dicts
    r_out_str = ''
    # query entries with same ndb_custom_key
    qry_entries_in_rand_dict = RandomDict.query(RandomDict.idx == ndb_custom_key, ancestor=ndb.Key(RandomDict, 'father'))  
    #print qry_entries_in_rand_dict.count()
    for qry in qry_entries_in_rand_dict:
        entty = qry.to_dict()
        #print entty['idx'] #entty['dict']
        rand_mapping_info_dict = entty['dict']
        rand_mapping_info_list = convert_shuffled_dict_to_str(rand_mapping_info_dict, c)
        rand_iter_numb = entty['entry_id'].split('~~~~')[-1] # last element from list
        rand_o_dir = rand_iter_numb + OUTPFILE
        # double check with if condition to see that the random dict are used for the correct run
        result = exec_core_microb_cmd(otu_table_biom, rand_o_dir, rand_mapping_info_list, c, group) 
        for r_frac_thresh , r_core_OTUs_biom in sorted(result['frac_thresh_core_OTUs_biom'].items(), key=lambda (key, value): int(key)): # return the items in sorted order
            r_OTUs , r_biom = r_core_OTUs_biom
            r_out_str += "Frac Threshold %s:\n%s\n%s\n\n" % (r_frac_thresh, ''.join(r_OTUs), r_biom)
            
            ndb_custom_key_r_frac_thres = ndb_custom_key + '~~~~' + r_frac_thresh
            ndb_custom_key_r_frac_thres_entry = ndb_custom_key_r_frac_thres + '~~~~' + rand_iter_numb


            Result_RandomDict(parent=ndb.Key(Result_RandomDict, 'fatherresults'), \
                     idx= ndb_custom_key, frac_thresh = ndb_custom_key_r_frac_thres, \
                     entry_id = ndb_custom_key_r_frac_thres_entry, \
                     otus = r_OTUs, biom = r_biom).put()


    '''
    ## using mapreduce to parallelize the core microbiome on random dicts
    shuffled_core_mic = ShuffleDictPipeline() 
    shuffled_core_mic.start()
    '''

 
    '''
    the following section compiles results from the Result Datatstore and calculates stats.
    '''
    sign_results = ''
    # compile results; print the number of random occurances for each true core microbiome otu (checks significance)
    for frac_s in [75, 80, 85, 90, 95, 100]:
        
        sign_results += 'Frac thresh %s\n' % str(frac_s)
        
        ndb_custom_key_qury_id = ndb_custom_key + '~~~~' + str(frac_s)
       
        qry_entries_in_result_rand_dict = Result_RandomDict.query(Result_RandomDict.frac_thresh == ndb_custom_key_qury_id, ancestor=ndb.Key(Result_RandomDict, 'fatherresults'))  
        #print qry_entries_in_result_rand_dict.count()
        # compile the results from randomization
        # this returns a list of list i.e. collects the unique set of core taxa from each randomized data
        taxons_ = list()
        for q in qry_entries_in_result_rand_dict:
            #print res
            res = q.to_dict()
            #print res['otus']
            taxons_.append(compile_results(res['otus']))
            
        all_results_otu_list = [item for sublist in taxons_ for item in sublist] #  taxons_ -> sublist -> item

        # calculate freq of otu being a core microb from the randomizations
        randomized_otus = calc_freq_elem_list(all_results_otu_list) # dict of otus and freq of occurance from random data
        '''
        #print the compiled results from randomized runs of core microbiomes
        for k, v in randomized_otus.items():
        	sign_results += "%s\t%s\n" % (k , v)
        '''

        # check significance
        #signif_core_microb_otu_list = list()
        for o in true_result_frac_thresh_otus_dict[str(frac_s)]:
            if o in randomized_otus:
                freq = int(randomized_otus[o])
                if freq/float(NTIMES) < 0.05:
                    otu = '%s\t%s\n' % (o, freq)
                    sign_results += otu
                    #signif_core_microb_otu_list.append(otu)
        

    return sign_results #r_out_str #result['otu_counts']

def convert_shuffled_dict_to_str(DICT, categ):
    file_str_list = list()
    file_str = "%s\t%s\n" % ('#SampleID' , categ)
    file_str_list.append(file_str)

    for k, v in DICT.items():
        for i in v:
            f_str = "%s\t%s\n" % (i , k)
            file_str_list.append(f_str)

    return file_str_list




MAIN_PAGE_HTML = """\
<html>
  <body>
    <form action="/sign" enctype="multipart/form-data" method="post">
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
	<div>
		<input type="submit" value="Send">
	</div>
    </form>
  </body>
</html>
"""


class MainPage(webapp2.RequestHandler):
    def get(self):
        self.response.write(MAIN_PAGE_HTML)

class Guestbook(webapp2.RequestHandler):
    def post(self):
        global DELIM, group, NTIMES, OUTPFILE , FRAC_SAMPLES
        DELIM = '\t'
    	factor, group = self.request.get('group').split(':')
    	NTIMES = int(self.request.get('random'))
    	datetime = strftime("%a-%d-%b-%Y-%I:%M:%S-%p", localtime())
    	## add a random alpha numeric to avoid conflict with another (simultaneous) user request.
    	OUTPFILE = self.request.get('content') + datetime
    
    	FRAC_SAMPLES = 100
        p_val_adj = self.request.get('pvaladjmethod')

        otu_table_biom = self.request.get('datafile') #automatically reads file
    	group_info = self.request.get('groupfile') #automatically reads file
        
        # check if biom or txt
        # later on use validate biom
        
        header_key="taxonomy"
        infile_txt = main_converter(otu_table_biom , header_key)
        
        dump_content = calc_significance(otu_table_biom, factor, group, group_info.split('\n'), p_val_adj)


        #dump_content = cgi.escape(otu_table_biom)  +  "\n" + cgi.escape(group_info) 
        #dump_content = infile_txt
        self.response.write('<html><body>You wrote:<pre>')
        self.response.write(dump_content)
        self.response.write('</pre></body></html>')


app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/sign', Guestbook),
], debug=True)


