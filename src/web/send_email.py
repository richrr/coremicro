# Copyright 2016, 2017 Richard Rodrigues, Nyle Rodgers, Mark Williams,
# Virginia Tech
#
# This file is part of Coremic.
#
# Coremic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Coremic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Coremic. If not, see <http://www.gnu.org/licenses/>.
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
        logging.error(result.error_info)
        logging.error(result.error_message)
    else:
        logging.info('Email sent to %s with subject %s', to_email, subj)
