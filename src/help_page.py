import webapp2
import run_config


class HelpPage(webapp2.RequestHandler):
    def get(self):
        template = run_config.JINJA_ENVIRONMENT.get_template('help.html')
        self.response.write(template.render())
