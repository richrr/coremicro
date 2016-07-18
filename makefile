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
