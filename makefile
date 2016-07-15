GAE_PROJECT_NAME = coremic2

deploy:
	appcfg.py -A $(GAE_PROJECT_NAME) update .

devel: 
	dev_appserver.py .

install: cleanlib
	wget -P lib/ https://github.com/biocore/biom-format/archive/1.2.0.zip --no-check-certificate
	unzip -q -d lib/ lib/1.2.0
	rm lib/1.2.0
	mv lib/biom-format-1.2.0/python-code/biom lib/biom
	rm -rf lib/biom-format-1.2.0/
	pip2 install -r requirements.txt -t lib/

cleanlib:
	rm -rf lib/*

clean: cleanlib
	rm *.pyc
