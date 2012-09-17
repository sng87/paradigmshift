#!/usr/bin/env python
"""snvScore.py: 

Usage:
  snvScore.py [options] snvFile proteinFa

Options:
  -q            run quietly
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import getopt, os, re, sys

verbose = True

def usage(code = 0):
    print __doc__
    if code != None:
        sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % (exitstatus)
        sys.exit(10)
    log("... done\n")

def factorial(n): 
    if n < 2:
        return(1)
    else:
        return(reduce(lambda x, y: x*y, xrange(2, int(n)+1)))

def binCDF(s, p, n):
    x = 1.0 - p
    a = n - s
    b = s + 1.0
    c = a + b - 1.0
    prob = 0
    for j in xrange(a, c + 1):
        prob += factorial(c) / (factorial(j)*factorial(c-j)) \
                * x**j * (1 - x)**(c-j)
    return(prob)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        log("ERROR: incorrect number of arguments", die = True)
    
    snvFile = args[0]
    proteinFa = args[1]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## reading protein files
    proteinLength = {}
    f = open(proteinFa, "r")
    for line in f:
        if line.isspace():
            continue
        if line.startswith(">"):
            pline = re.split("\t", line.rstrip("\r\n"))
            gene = pline[0].lstrip(">")
            length = int(pline[1])
            proteinLength[gene] = length
    f.close()
    
    ## reading snvs
    snvMap = {}
    f = open(snvFile, "r")
    for line in f:
        if line.isspace():
            continue
        pline = re.split("\t", line.rstrip("\r\n"))
        gene = pline[0]
        position = int(pline[2][1:-1])
        if gene not in snvMap:
            snvMap[gene] = {}
            snvMap[gene]["total"] = 0
        if position not in snvMap[gene]:
            snvMap[gene][position] = 1
            snvMap[gene]["total"] += 1
        else:
            snvMap[gene][position] += 1
            snvMap[gene]["total"] += 1
    f.close()
    
    ## computing cross-protein cdf
    for gene in proteinLength.keys():
        if gene not in snvMap:
            continue
        f = open("%s.cdf.txt" % (gene), "w")
        x = 0.0
        y = 0.0
        hits = 0
        f.write("%s\t%s\n" % (x, y))
        for position in range(1, proteinLength[gene]+1):
            x = float(position)/(proteinLength[gene]+1)
            y = float(hits)/snvMap[gene]["total"]
            if position in snvMap[gene]:
                f.write("%s\t%s\n" % (x, y))
                hits += snvMap[gene][position]
                y = float(hits)/snvMap[gene]["total"]
                f.write("%s\t%s\n" % (x, y))
        f.write("1.0\t1.0\n")
        f.close()
        
    ## compute per-position scores
    windowWidth = 100
    for gene in proteinLength.keys():
        if gene not in snvMap:
            continue
        f = open("%s.score.txt" % (gene), "w")
        for origin in snvMap[gene].key():
            windowStart = origin - windowWidth
            if windowStart < 1:
                windowStart = 1
            windowStop = origin + windowWidth
            if windowStop > proteinLength[gene]:
                windowStop = proteinLength[gene]
            windowSize = windowStop - windowStart
            nMut = 0
            for position in range(windowStart, windowStop+1):
                if position in snvMap[gene]:
                    nMut += snvMap[gene][position]
            pScore = 1 - binCDF(nMut - 1, float(windowSize)/float(proteinLength[gene]), snvMap[gene]["total"])
            f.write("%s\t%s\t%s\n" % (origin, snvMap[gene][position], pScore))
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
