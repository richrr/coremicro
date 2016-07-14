import os
import datetime


# Whether this is running in production or development
IS_PRODUCTION = os.getenv('SERVER_SOFTWARE',
                          '').startswith('Google App Engine/')

# How many parallel pipes to have processing the randomized data
MAX_NUM_PARALLEL = 8

# If it looks like the task will run longer than this to do another itteration
# start a new task
MAX_RUNNING_TIME = datetime.timedelta(minutes=9)

# The fractions to run the analysis at
FRACS = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]

# The delimiter used in the groupfile
DELIM = '\t'
