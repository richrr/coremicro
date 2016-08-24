#!/usr/bin/env python
import webapp2

from main_page import MainPage

app = webapp2.WSGIApplication([
    ('/', MainPage)
], debug=True)
