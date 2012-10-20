#!/usr/bin/env python
"""mutSummary.py:

Usage:
mutSummary.py [options] maf-file

Options:
   -q          run quietly, don't output status
"""
# Written by: Sam Ng
# Last updated: 3/11/11
import getopt, os, os.path, random, re, sys, time, glob
from copy import deepcopy

from mSHIFT import *
import mData, mCalculate, mPathway

verbose = True

scriptDir = os.path.realpath(os.path.dirname(sys.argv[0]))
scriptDir += "/"

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if (verbose):
        sys.stderr.write(msg)
    if (die):
        sys.exit(1)

def syscmd(cmd):
    log("\t" + cmd + "\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "ERROR: Failed with exit status %i" % exitstatus
        sys.exit(10)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    if len(args) != 1:
        print "incorrect number of arguments"
        usage(1)
    
    mutFile = args[0]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## read mafFile
    (mutMap, mutClass, allSamples) = rMAF(mutFile)
    
    ## write mutationSummary
    mutGenes = mutMap.keys()
    mutGenes.sort(lambda x,y: cmp(len(mutMap[y]), len(mutMap[x])))
    f = open("mutationSummary.tab", "w")
    for i in mutGenes:
        mutSamples = list()
        for j in mutMap[i]:
            mutSamples.append(j[0:16])
        f.write("%s\t%s\t%s\n" % (i, len(mutMap[i]), ",".join(mutSamples)))
    f.close()
    mData.wList("include.samples", allSamples)
    
    ## write mutation.tab
    f = open("mutation.tab", "w")
    f.write("id\t%s\n" % ("\t".join(allSamples)))
    for i in mutGenes:
        f.write("%s" % (i))
        for j in allSamples:
            if j in mutMap[i]:
                f.write("\t1")
            else:
                f.write("\t0")
        f.write("\n")
    f.close()
    
    ## temporarily hardcoded output of include.samples for TP53
    if "TP53" in mutClass:
        mData.wList("include.samples.TP53.truncating" , list(set(allSamples) - set(mutClass["TP53"]["missense"]) - set(mutClass["TP53"]["silent"]) - set(mutClass["TP53"]["other"])))
        mData.wList("include.samples.TP53.missense" , list(set(allSamples) - set(mutClass["TP53"]["truncating"]) - set(mutClass["TP53"]["silent"]) - set(mutClass["TP53"]["other"])))
        
    
if __name__ == "__main__":
    main(sys.argv[1:])
