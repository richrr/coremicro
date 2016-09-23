import logging
import mailjet_rest
import requests_toolbelt.adapters.appengine

from email_config import API_KEY, SECRET_KEY, FROM_EMAIL

# Use the App Engine requests adapter to allow the requests library to be
# used on App Engine.
requests_toolbelt.adapters.appengine.monkeypatch()


def send_email(subj, msg, to_email, attachments=[]):
    mailjet = mailjet_rest.Client(auth=(API_KEY, SECRET_KEY))
    data = {
        'FromEmail': FROM_EMAIL,
        'FromName': 'Coremic',
        'Subject': subj,
        'Text-part': msg,
        'Recipients': [
            {
                'Email': to_email,
            }
        ],
        'Attachments': attachments,
    }
    result = mailjet.send.create(data=data)
    if result.status_code != 200:
        logging.error('Unable to sucesfully email results\n' +
                      'Status code of %s' % result.status_code)
