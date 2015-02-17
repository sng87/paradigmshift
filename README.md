PARADIGM-SHIFT: Prediction of Functional Impact Using Pathway Integration
========

Current Version 
--------

1.0

Authors
--------

Sam Ng and Joshua M. Stuart


Requirements
--------

* [python](http://www.python.org/) >= 2.7
   * [scipy](http://www.scipy.org/) >= 0.12.0
   * [numpy](http://numpy.scipy.org/)
   * [pandas](http://pandas.pydata.org/)
   * [networkx](http://networkx.github.io/)

Installation
-------

- Install dependencies
- Download the paradigmshift repository to the desired location
- Run "make" in paradigmshift/ to generate source files; may also require downloading the paradigm-scripts and pathmark-scripts repository (see https://github.com/ucscCancer/paradigm-scripts and https://github.com/ucscCancer/pathmark-scripts)
- Source init files for paradigmshift, pathmark-scripts and paradigm-scripts (init.sh for bash and init.csh for csh)
- Run code on example data in examples/ with "make"

Command-Line
------
```
paradigmSHIFT.py [options] paradigm_directory analysis_file

paradigm_directory - directory or zip of the paradigm cohort to perform analysis on
analysis_file - file containing events to run paradigm-shift analysis on

-w work_directory - directory to perform work in
-c config_file - file with user defined arguments
-s include_samples - one column file of samples to analyze
-f include_features - one column file of features to analyze
-n null_size - number of null samples to estimate shift significance
-b batch_size - number of samples in each paradigm run per compute node
-p pathway_file - path to paradigm pathway file
-y - this flag must be enabled when using the publicly available version of paradigm
-z seed - fix random seed
```

Folders
------
* bin : executables
* examples : LUSC and LUAD inputs for demonstration purposes

Contact
------
Feature requests, comments and requests for clarification should all be sent to the author at <sng@soe.ucsc.edu>. 
I will try to respond quickly to all requests, so feel free to email me!
