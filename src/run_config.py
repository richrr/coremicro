import os

# Whether this is running in production or development
IS_PRODUCTION = os.getenv('SERVER_SOFTWARE',
                          '').startswith('Google App Engine/')

# The fractions to run the analysis at
FRACS = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]

# The delimiter used in the groupfile
DELIM = '\t'
