# About
Coremic is a tool finding the core microbiome of some group from user suplied input data. It runs on Google App Engine.

# Dependency Instalation
After cloning, the dependencies can be installed by running `make install` from the main directory of the project. This assumes that you have make, wget, unzip, and pip installed. Alternatively one can also install the dependencies manually. To do this first download [biom-format 1.2.0](https://github.com/biocore/biom-format/archive/1.2.0.zip) and unzip it. Then place the biom/ directory (found at biom-format-1.2.0/python-code/biom inside the unziped directory) in the lib/ directory of the project. Finally run `pip2 install -t lib/ -r requirements.txt` from the main directory of the project to install the remaining dependencies.

# Testing and Deploying
`make devel` will start the development server, and `make deploy` will deploy to GAE. You must have first installed all dependencies for these to work. Before running `make deploy` change `GAE_PROJECT_NAME` in the makefile to your GAE project name.

# Email Configuration
The web version of this tool delivers its results via email. To configure this create `src/web/email_config.py`, which should look like the `example_email_config.py` in that same directory, replacing the keys and from address with the keys for your Mailjet account and the email with the email you wish to send the results from (make sure this is set up with Mailjet as well)!