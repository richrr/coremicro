import webapp2
from main_page import MainPage
from help_page import HelpPage

# Set up web server
app = webapp2.WSGIApplication([
    ('/', MainPage),
    ('/help', HelpPage),
], debug=True)
