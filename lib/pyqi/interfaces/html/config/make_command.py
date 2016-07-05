#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__credits__ = ["Evan Bolyen"]

from pyqi.core.interfaces.html import (HTMLInputOption, HTMLDownload, HTMLPage)
from pyqi.core.interfaces.optparse.input_handler import string_list_handler
from pyqi.core.interfaces.html.output_handler import newline_list_of_strings
from pyqi.core.command import (make_command_in_collection_lookup_f,
    make_command_out_collection_lookup_f)
from pyqi.commands.make_command import CommandConstructor

cmd_in_lookup = make_command_in_collection_lookup_f(CommandConstructor)
cmd_out_lookup = make_command_out_collection_lookup_f(CommandConstructor)

inputs = [
    HTMLInputOption(Parameter=cmd_in_lookup('name')),
    HTMLInputOption(Parameter=cmd_in_lookup('author')),
    HTMLInputOption(Parameter=cmd_in_lookup('email')),
    HTMLInputOption(Parameter=cmd_in_lookup('license')),
    HTMLInputOption(Parameter=cmd_in_lookup('copyright')),
    HTMLInputOption(Parameter=cmd_in_lookup('version'), Name='command-version'),
    HTMLInputOption(Parameter=cmd_in_lookup('credits'),
                   Handler=string_list_handler,
                   Help='comma-separated list of other authors'),
    HTMLInputOption(Parameter=cmd_in_lookup('test_code'),
                   Type=bool,
                   Help='Should a stubbed out python test file be generated instead'),
    HTMLInputOption(Parameter=None,
                   Name='download-file',
                   Required=True,
                   Help='The name of the file to download which contains generated Python code. (e.g. my_command)')
]

outputs = [ HTMLDownload(Parameter=cmd_out_lookup('result'),
                   Handler=newline_list_of_strings,
                   FilenameLookup='download-file',
                   FileExtension='.py') ]

#Comment out the above and uncomment the below for an example of a page.

#    [ HTMLPage(Parameter=cmd_out_lookup('result'),
#              Handler=newline_list_of_strings) ]
