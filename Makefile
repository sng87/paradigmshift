THISDIR = ${CURDIR}
THISOS = ${shell uname -s}

PARADIGM_GIT = https://github.com/ucscCancer/paradigm-scripts.git

all : init.sh init.csh

init.sh :
	echo export PATH=${THISDIR}/bin:\$${PATH} > init.sh
	echo if [ -n "\$${PYTHONPATH+x}" ] >> init.sh
	echo then >> init.sh
	echo export PYTHONPATH=${THISDIR}/bin:\$${PYTHONPATH} >> init.sh
	echo else >> init.sh
	echo export PYTHONPATH=${THISDIR}/bin >> init.sh
	echo fi >> init.sh

init.csh :
	echo setenv PATH ${THISDIR}/bin:\$${PATH} > init.csh
	echo if \$$?PYTHONPATH then >> init.csh
	echo setenv PYTHONPATH ${THISDIR}/bin:\$${PYTHONPATH} >> init.csh
	echo else >> init.csh
	echo setenv PYTHONPATH ${THISDIR}/bin >> init.csh
	echo endif >> init.csh

paradigm-scripts :
	if [ ! -d '../paradigm-scripts' ]; then \
		cd ..; git clone ${PARADIGM_GIT}; cd paradigm-scripts; make; \
	fi
	ln -s ../paradigm-scripts paradigm-scripts

clean :
	rm -f init.sh init.csh
	if [ -d 'example' ]; then \
		cd example; make clean; \
	fi
