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

verbose = True

scriptDir = os.path.realpath(os.path.dirname(sys.argv[0]))
scriptDir += "/"

truncList = ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"]
missList = ["Missense_Mutation"]


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

def wList(outf, outList):
     """write 1 column list"""
     f = open(outf, "w")
     for i in outList:
         f.write("%s\n" % (i))
     f.close()

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
    wList("include.samples", allSamples)
    
    ## write mutation.tab
    f = open("mutation.tab", "w")
    f.write("id\t%s\n" % ("\t".join(allSamples)))
    for i in mutGenes:
        f.write("%s" % (i))
        for j in allSamples:
            types = []
            for type in mutClass[i].keys():
                if j in mutClass[i][type]:
                    types.append(j)
            if len(types) > 0:
                f.write("\t%s" % (",".join(types)))
            else:
                f.write("\tNone")
        f.write("\n")
    f.close()
    
    ## temporarily hardcoded output of include.samples
    focusGene = "TP53"
    if focusGene in mutClass:
        truncatingSamples = []
        for type in truncList:
            if type in mutClass[focusGene]:
                for sample in mutClass[focusGene][type]:
                    if sample not in truncatingSamples:
                        truncatingSamples.append(sample)
        missenseSamples = []
        for type in missList:
            if type in mutClass[focusGene]:
                for sample in mutClass[focusGene][type]:
                    if sample not in missenseSamples:
                        missenseSamples.append(sample)    
        wList("%s.truncating.tab" % (focusGene), truncatingSamples)
        wList("%s.missense.tab" % (focusGene), missenseSamples)
        
    
if __name__ == "__main__":
    main(sys.argv[1:])
