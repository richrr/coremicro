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
GAE_PROJECT_NAME = coremic2
BIOM_VERSION = 1.1.2

deploy:
	appcfg.py -A $(GAE_PROJECT_NAME) update .

devel: 
	dev_appserver.py .

install: cleanlib
	wget -P lib/ https://github.com/biocore/biom-format/archive/$(BIOM_VERSION).zip --no-check-certificate
	unzip -q -d lib/ lib/$(BIOM_VERSION)
	rm lib/$(BIOM_VERSION)
	mv lib/biom-format-$(BIOM_VERSION)/python-code/biom lib/biom
	rm -rf lib/biom-format-$(BIOM_VERSION)/
	pip2 install -r requirements.txt -t lib/

cleanlib:
	rm -rf lib/*

clean: cleanlib
	rm *.pyc
