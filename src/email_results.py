from google.appengine.api import mail
from google.appengine.api.app_identity import get_application_id

import logging


def send_results_as_email(params, attachments):
    subj = "Your data from %s with name %s has been processed" % (
        params['timestamp'], params['run_name'])
    msg_str = """
Dear User:

Your data has been processed and is attached. Thanks for using this tool.

Please email us if you have any questions.

The Core Microbiome Team

"""
    msg_str += params['user_args']
    message = make_base_email(subj, params['to_email'], msg_str)
    message.body = msg_str
    message.attachments = attachments

    message.send()


def send_error_as_email(params, error):
    logging.warn(error)
    subj = 'There was an error in processing your data from %s with name %s'\
           % (params['timestamp'], params['run_name'])
    msg_str = """
Dear User:

There was an error in processing your data. The error is listed below.

Please email us if you have any questions.

The Core Microbiome Team

"""
    msg_str += params['user_args']
    msg_str += '\n' + error
    message = make_base_email(subj, params['to_email'], msg_str)
    message.send()


def make_base_email(subj, to_email, msg):
    app_id = get_application_id()
    return mail.EmailMessage(
        sender="Core Microbiome <do-not-reply@%s.appspotmail.com>" % app_id,
        subject=subj,
        to='User <' + to_email + '>',
        body=msg)
