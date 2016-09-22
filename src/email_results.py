from time import strptime, mktime
from datetime import datetime
import logging
import mailjet_rest
import requests_toolbelt.adapters.appengine

from email_config import *

# Use the App Engine requests adapter to allow the requests library to be
# used on App Engine.
requests_toolbelt.adapters.appengine.monkeypatch()


def send_results_as_email(params, attachments):
    subj = "Your data with name %s has been processed" % params['run_name']
    msg_str = """
Dear User:

Your data has been processed and is attached. Thanks for using this tool.

Please email us if you have any questions.

The Core Microbiome Team

"""
    msg_str += params['user_args']
    elapsed_time = datetime.now() - datetime.fromtimestamp(mktime(
        strptime(params['timestamp'], '%a-%d-%b-%Y-%I:%M:%S-%p')))
    logging.info('Elapsed time: %d.%06d seconds' % (elapsed_time.seconds,
                                                    elapsed_time.microseconds))
    msg_str += 'Elapsed time: %d.%06d seconds\n' % (elapsed_time.seconds,
                                                    elapsed_time.microseconds)

    mailjet = mailjet_rest.Client(auth=(API_KEY, SECRET_KEY))
    data = {
        'FromEmail': FROM_EMAIL,
        'FromName': 'Coremic',
        'Subject': subj,
        'Text-part': msg_str,
        'Recipients': [
            {
                'Email': params['to_email'],
            }
        ],
        'Attachments': attachments,
    }
    result = mailjet.send.create(data=data)
    if result.status_code != 200:
        logging.error('Unable to sucesfully email results\n' +
                      'Status code of %s' % result.status_code)


def send_error_as_email(params, error):
    logging.warn(error)
    subj = 'There was an error in processing your data with name %s'\
           % params['run_name']
    msg_str = '''
Dear User:

There was an error in processing your data. The error is listed below.

Please email us if you have any questions.

The Core Microbiome Team

'''
    msg_str += params['user_args']
    elapsed_time = datetime.now() - datetime.fromtimestamp(mktime(
        strptime(params['timestamp'], '%a-%d-%b-%Y-%I:%M:%S-%p')))
    logging.info('Elapsed time: %d.%06d seconds' % (elapsed_time.seconds,
                                                    elapsed_time.microseconds))
    msg_str += 'Elapsed time: %d.%F6d seconds\n' % (elapsed_time.seconds,
                                                    elapsed_time.microseconds)
    msg_str += '\n' + error

    mailjet = mailjet_rest.Client(auth=(API_KEY, SECRET_KEY))
    data = {
        'FromEmail': 'donotreply@coremic2.appspot.com',
        'FromName': 'Coremic',
        'Subject': subj,
        'Text-part': msg_str,
        'Recipients': [
            {
                'Email': params['to_email'],
            }
        ],
    }
    result = mailjet.send.create(data=data)
    if result.status_code != 200:
        logging.error('Unable to sucesfully email results\n' +
                      'Status code of %s' % result.status_code)
