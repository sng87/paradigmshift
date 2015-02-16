THISDIR = ${CURDIR}

PARADIGM_GIT = git://github.com/ucscCancer/paradigm-scripts.git
PATHMARK_GIT = git://github.com/ucscCancer/pathmark-scripts.git

all : init.sh init.csh

init.sh : bin/circlePlot.py
	echo export PATH=${THISDIR}/bin:\$${PATH} > init.sh
	echo if [ -n "\$${PYTHONPATH+x}" ] >> init.sh
	echo then >> init.sh
	echo export PYTHONPATH=${THISDIR}/bin:\$${PYTHONPATH} >> init.sh
	echo else >> init.sh
	echo export PYTHONPATH=${THISDIR}/bin >> init.sh
	echo fi >> init.sh

init.csh : bin/circlePlot.py
	echo setenv PATH ${THISDIR}/bin:\$${PATH} > init.csh
	echo if \$$?PYTHONPATH then >> init.csh
	echo setenv PYTHONPATH ${THISDIR}/bin:\$${PYTHONPATH} >> init.csh
	echo else >> init.csh
	echo setenv PYTHONPATH ${THISDIR}/bin >> init.csh
	echo endif >> init.csh

bin/circlePlot.py : paradigm-scripts
	ln -s $$PWD/paradigm-scripts/bin/circlePlot.py bin/circlePlot.py

paradigm-scripts :
	if [ ! -d '../paradigm-scripts' ]; then \
		cd ..; git clone ${PARADIGM_GIT}; cd paradigm-scripts; make; \
	fi
	ln -s ../paradigm-scripts paradigm-scripts

pathmark-scripts :
	if [ ! -d '../pathmark-scripts' ]; then \
		cd ..; git clone ${PATHMARK_GIT}; cd pathmark-scripts; make; \
	fi
	ln -s ../pathmark-scripts pathmark-scripts

clean :
	rm -f init.sh init.csh bin/circlePlot.py paradigm-scripts pathmark-scripts
	if [ -d 'example' ]; then \
		cd example; make clean; \
	fi
