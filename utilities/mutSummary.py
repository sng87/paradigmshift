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

def rMAF(inf, delim = "\t"):
    """read .maf format file"""
    mutData = {}
    mutClass = {}
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    pline = line.rstrip().split(delim)
    hugoCol = -1
    tumorCol = -1
    classCol = -1
    for i, j in enumerate(pline):
        if j == "Hugo_Symbol":
            hugoCol = i
        elif j == "Tumor_Sample_Barcode":
            tumorCol = i
        elif j == "Variant_Classification":
            classCol = i
    samples = []
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split(delim)
        if pline[hugoCol] not in mutData:
            mutData[pline[hugoCol]] = []
            mutClass[pline[hugoCol]] = {}
        mutData[pline[hugoCol]].append(pline[tumorCol])
        if classCol != -1:
            if pline[classCol] not in mutClass[pline[hugoCol]]:
                mutClass[pline[hugoCol]][pline[classCol]] = []
            mutClass[pline[hugoCol]][pline[classCol]].append(pline[tumorCol])
        if pline[tumorCol] not in samples:
            samples.append(pline[tumorCol])
    f.close()
    return(mutData, mutClass, samples)

def wList(outf, outList):
     """write 1 column list"""
     f = open(outf, "w")
     for i in outList:
         f.write("%s\n" % (i))
     f.close()

def formatSample(sample):
    """this function may need to be modified depending on the project"""
    return(sample[0:15])

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
        mutSamples = []
        for j in mutMap[i]:
            mutSamples.append(formatSample(j))
        f.write("%s\t%s\t%s\n" % (i, len(mutMap[i]), ",".join(mutSamples)))
    f.close()
    wList("include.samples", [formatSample(i) for i in allSamples])
    
    ## write mutation.list_t
    f = open("mutation.list_t", "w")
    for i in mutGenes:
        mutSamples = []
        for j in mutMap[i]:
            mutSamples.append(formatSample(j))
        f.write("%s\t%s\n" % (i, "\t".join(mutSamples)))
    f.close()
    
    ## write mutation.tab
    f = open("mutation.tab", "w")
    f.write("id\t%s\n" % ("\t".join(allSamples)))
    for i in mutGenes:
        f.write("%s" % (i))
        for j in allSamples:
            types = []
            for type in mutClass[i].keys():
                if j in mutClass[i][type]:
                    types.append(type)
            if len(types) > 0:
                f.write("\t%s" % (",".join(types)))
            else:
                f.write("\tNone")
        f.write("\n")
    f.close()
    
    ## write classifications
    f = open("mutationClass.tab", "w")
    for gene in mutClass:
        truncSamples = []
        missSamples = []
        negSamples = [formatSample(i) for i in allSamples]
        for type in truncList:
            if type in mutClass[gene]:
                for sample in mutClass[gene][type]:
                    if formatSample(sample) not in truncSamples:
                        truncSamples.append(formatSample(sample))
        for type in missList:
            if type in mutClass[gene]:
                for sample in mutClass[gene][type]:
                    if formatSample(sample) not in missSamples:
                        missSamples.append(formatSample(sample))
        negSamples = list(set(negSamples) - (set(truncSamples) | set(missSamples)))
        f.write("%s\t%s\t%s\t%s\t%s\n" % (gene, len(truncSamples), ",".join(truncSamples), ",".join(missSamples), ",".join(negSamples)))
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
