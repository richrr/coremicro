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
import webapp2

from process_data import ProcessData
from process_results import ProcessResults
from main_page import MainPage

# to do the parallelism: async, futures or mapreduce
# https://cloud.google.com/appengine/docs/python/datastore/async#Working_with_the_Async_Datastore_API
# futures is for python 3, so cannot be used for my 2.7
# mapreduce should work; if it doesn't, try the async datastore api


# http://stackoverflow.com/questions/11849456/how-to-filter-datastore-data-before-mapping-to-cloud-storage-using-the-mapreduce
# http://stackoverflow.com/questions/23508116/appengine-mapreduce-how-to-filter-structuredproperty-while-using-datastore-input


app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/process_data', ProcessData),
    ('/process_results', ProcessResults),
], debug=True)
