# Copyright 2016 Richard Rodrigues, Nyle Rodgers, Mark Williams, Virginia Tech
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
import os
import jinja2

# Whether this is running in production or development
IS_PRODUCTION = os.getenv('SERVER_SOFTWARE',
                          '').startswith('Google App Engine/')

# The fractions to run the analysis at
FRACS = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]

# The delimiter used in the groupfile
DELIM = '\t'

# For templates
JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__),
                                                '../templates/')),
    extensions=['jinja2.ext.autoescape'],
    autoescape=True)
