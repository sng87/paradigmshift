#!/usr/bin/env python
"""median-center.py: performs gene-wise normalization of a tab-file (vs. normal if -n is specified)

Usage:
  median-center.py [options] input-tab output-tab

Options:
  -m            use mean instead
  -n str        normal-tab file of normal samples
  -s str        use samples in list only
  -q            run quietly
"""
## Written By: Sam Ng
## Last Updated: 8/12/11
import os, os.path, sys, getopt, re

useMean = False

verbose = True

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log("... done\n")

def rList(inf, header = False):
    """read 1 column list"""
    inList = []
    f = open(inf, "r")
    if header:
        f.readline()
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\t\r\n")
        inList.append(line)
    f.close()
    return(inList)

def rCRSData(inf, delim = "\t", retFeatures = False):
    """reads .tsv into a [col][row] dictionary"""
    inData = dict()
    colFeatures = []
    rowFeatures = []
    ## read header
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    colFeatures = pline[1:]
    for i in colFeatures:
        if i not in inData:
            inData[i] = dict()
    ## read data
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        rowFeatures.append(pline[0])
        if len(pline) != (1+len(colFeatures)):
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        for i in range(len(colFeatures)):
            if pline[i+1] == "":
                inData[colFeatures[i]][pline[0]] = "NA"
            else:            
                inData[colFeatures[i]][pline[0]] = pline[i+1]
    f.close()
    if retFeatures:
        return(inData, colFeatures, rowFeatures)
    else:
        return(inData)

def wCRSData(outf, outData, delim = "\t", useCols = None, useRows = None):
    """write [col][row] dictionary to .tsv"""
    ## get colFeatures and rowFeatures
    if useCols is None:
        colFeatures = outData.keys()
    else:
        colFeatures = useCols
    if useRows is None:
        rowFeatures = outData[colFeatures[0]].keys()
    else:
        rowFeatures = useRows
    ## write header
    f = open(outf, "w")
    f.write("id")
    for i in colFeatures:
        f.write("\t%s" % (i))
    f.write("\n")
    for i in rowFeatures:
        f.write("%s" % (i))
        for j in colFeatures:
            if j in outData:
                if i in outData[j]:
                    f.write("\t%s" % (outData[j][i]))
                else:
                    f.write("\tNA") 
            else:
                f.write("\tNA")
        f.write("\n")
    f.close()

def floatList(inList):
    """returns only numeric elements of a list"""
    outList = []
    for i in inList:
        try:
            fval = float(i)
            if fval != fval:
                raise ValueError
            outList.append(fval)
        except ValueError:
            continue
    return(outList)

def median(inList):
    """calculates median"""
    cList = floatList(inList)
    cList.sort()
    if len(cList) == 0:
        median = "NA"
    else:
        if len(cList)%2 == 1:
            median = cList[len(cList)/2]
        else:
            median = (cList[len(cList)/2]+cList[(len(cList)/2)-1])/2.0
    return(median)

def mean(inList, null = "NA"):
    """Calculates mean"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = null
    else:
        mean = sum(fList)/len(fList)
    return (mean)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "mn:s:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        print "incorrect number of arguments"
        usage(1)
    
    inf = args[0]
    outf = args[1]
    
    normFile = None
    sampleFile = None
    includeSamples = []
    global verbose, useMean
    for o, a in opts:
        if o == "-m":
            useMean = True
        elif o == "-n":
            normFile = a
        elif o == "-s":
            sampleFile = a
        elif o == "-q":
            verbose = False
    
    ## read files
    (sampleData, cols, rows) = rCRSData(inf, retFeatures = True)
    if normFile is not None:
        normalData = rCRSData(normFile)
    if sampleFile is not None:
        includeSamples = rList(sampleFile)
 
    ## compute-medians
    medianMap = {}
    for i in rows:
        vals = []
        if normFile is not None:
            for j in normalData.keys():
                if i not in normalData[j]:
                    log("ERROR: Missing feature %s in normal matrix\n" % (i), die = True)
                if sampleFile is None:
                    vals.append(normalData[j][i])
                elif j in includeSamples:
                    vals.append(normalData[j][i])
        else:
            for j in cols:
                if sampleFile is None:
                    vals.append(sampleData[j][i])
                elif j in includeSamples:
                    vals.append(sampleData[j][i])
        if useMean:
            medianMap[i] = mean(vals)
        else:
            medianMap[i] = median(vals)
    
    ## median-center data
    outData = {}
    for i in cols:
        outData[i] = {}
        for j in rows:
            try:
                outData[i][j] = float(sampleData[i][j]) - medianMap[j]
            except ValueError:
                outData[i][j] =  "NA"
    wCRSData(outf, outData, useCols = cols, useRows = rows)

if __name__ == "__main__":
    main(sys.argv[1:])
