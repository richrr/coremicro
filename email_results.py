from google.appengine.api import mail


# the user id needs to be changed to that input by the user
def send_results_as_email(timestmp, user_args, msg, to_email):
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
    message.attachments = ['results.tsv', msg.encode('utf-8')]
    message.send()