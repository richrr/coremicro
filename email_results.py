from google.appengine.api import mail

import logging


# the user id needs to be changed to that input by the user
def send_results_as_email(timestmp, user_args, results, tree,
                          to_email,):
    subj = "Your data from %s has been processed" % timestmp
    message = mail.EmailMessage(
        sender="Core microbiome Support <richieangel@gmail.com>",
        subject=subj)

    # "User <richrr@vt.edu>" #Albert Johnson <Albert.Johnson@example.com>
    email_to = 'User <' + to_email + '>'
    message.to = email_to

    msg_str = """
Dear User:

Your data has been processed and is attached. Thanks for using this tool.

Please email us if you have any questions.

The Core Microbiome Team

"""
    msg_str += user_args
    message.body = msg_str
    message.attachments = [
        ('results.tsv', results.encode('utf-8')),
        ('results.nh', tree.encode('utf-8'))
    ]
    message.send()


def send_error_as_email(timestmp, user_args, error, to_email):
    logging.warn(error)
    subj = 'There was an error in processing your data from %s'\
           % timestmp
    message = mail.EmailMessage(
        sender="Core microbiome Support <richieangel@gmail.com>",
        subject=subj)

    # "User <richrr@vt.edu>" #Albert Johnson <Albert.Johnson@example.com>
    email_to = 'User <' + to_email + '>'
    message.to = email_to

    msg_str = """
Dear User:

There was an error in processing your data. The error is listed below.

Please email us if you have any questions.

The Core Microbiome Team

"""
    msg_str += user_args
    msg_str += '\n' + error
    message.body = msg_str
    message.send()
