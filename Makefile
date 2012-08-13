THISDIR = $(CURDIR)
THISOS = $(shell uname -s)

init.sh :
	echo \
	export PATH=$(THISDIR)/bin:$(THISDIR)/utilities:\$${PATH} > init.sh
	echo \
	if [ -n "\$${PYTHONPATH+x}" ] >> init.sh
	echo \
	then >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	else >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR):$(THISDIR)/bin >> init.sh
	echo \
	fi >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:$(THISDIR)/utilities:\$${PATH} > init.csh
	echo \
	if \$$?PYTHONPATH then >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR):$(THISDIR)/bin:\$${PYTHONPATH} >> init.csh
	echo \
	else >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR):$(THISDIR)/bin >> init.csh
	echo \
	endif >> init.csh

../paradigm-scripts :
	cd ..; git clone git://github.com/sng87/paradigm-scripts.git
	cd ../paradigm-scripts; make

clean :
	rm -f init.sh init.csh
	cd test; make clean
