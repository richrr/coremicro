#!/usr/bin/env python
import webapp2

from process_data import ProcessData
from main_page import MainPage

app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/process_data', ProcessData),
], debug=True)
