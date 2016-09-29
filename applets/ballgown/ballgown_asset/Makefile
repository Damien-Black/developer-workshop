all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo apt-get --yes update
	sudo apt-get remove --auto-remove --yes r-base
	sudo apt-get install --yes --force-yes r-base
	sudo apt-get install --yes --force-yes r-base-dev
	R -e 'source("http://bioconductor.org/biocLite.R"); biocLite(); biocLite("ballgown")'
