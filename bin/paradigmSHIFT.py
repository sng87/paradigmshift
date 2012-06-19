#!/usr/bin/env python
"""paradigmSHIFT.py: a script for identifying functionally important mutations across cancer cohorts

Usage:
  paradigmSHIFT.py [options]

Optional Configuration Files:
  
"""
## Written By: Sam Ng
import math, os, sys, random, re
from copy import deepcopy

import mData, mCalculate
from mSHIFT import *

from optparse import OptionParser
from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

randSeed = random.randint(0,9999999999999999)
f = open("seed.log", "w")
f.write("%s\n" % (randSeed))
f.close()

## default variables
nNulls = 10
mFolds = 5
rRepeats = 1
alpha = 0.05
usePerm = False
useParadigm = False

## params variables
mutThresh = 5
signalMethod = 0
cohortName = "default"
distanceParams = [2]
tstatParams = [1.0, 2.0]
incrParams = [0.1]
methodParams = ["vsMax"]
if os.path.exists("mut.cfg"):
    f = open("mut.cfg", "r")
    x = f.readlines()
    for line in x:
        if line.isspace():
            continue
        pline = re.split("\t", line.rstrip("\n"))
        if line.startswith("mutThresh"):
            mutThresh = int(pline[1])
        elif line.startswith("signalMethod"):
            signalMethod = int(pline[1])
        elif line.startswith("cohortName"):
            cohortName = pline[1]
        elif line.startswith("distanceParams"):
            distanceParams = [int(i) for i in re.split(",", pline[1])]
        elif line.startswith("tstatParams"):
            tstatParams = [float(i) for i in re.split(",", pline[1])]
        elif line.startswith("incrParams"):
            incrParams = [float(i) for i in re.split(",", pline[1])]
        elif line.startswith("methodParams"):
            methodParams = [i for i in re.split(",", pline[1])]
        elif line.startswith("vsPerm"):
            usePerm = True
        elif line.startswith("vsNorm"):
            usePerm = False
    f.close()

## executables and directories
paradigmExec = "/hive/users/sng/bin/paradigm-scripts/exe/paradigm"
mswarmExec = "/hive/users/sng/bin/paradigm-scripts/bin/mergeSwarmFiles.py"
mmergeExec = "/hive/users/sng/bin/paradigm-scripts/bin/merge_merged.py"
circleExec = "circlePlot.py"
medianExec = "median-center.py"
cytowebExec = "layout-cytoscapeweb.py"
htmlDir = "/hive/users/sng/.html"

## check files
pathFile = None
cnvFile = None
expFile = None
assert os.path.exists("paradigm")
assert os.path.exists("paradigm/clusterFiles")
for i in os.listdir("paradigm/clusterFiles"):
    if i.endswith("pathway.tab"):
        pathFile = i
    elif i.endswith("CNV.tab"):
        cnvFile = i
    elif (i.endswith("Expression.tab") | i.endswith("Expression.vCohort.tab") | i.endswith("Expression.vNormal.tab")):
        expFile = i
assert (pathFile != None)
assert (cnvFile != None)
assert (expFile != None)

## read global pathway
(gNodes, gInteractions) = rPathway("paradigm/clusterFiles/%s" % (pathFile))
gfPathway = flattenPathway(Pathway(gNodes, gInteractions))
#gPathway = Pathway(gNodes, gInteractions)
gPathway = gfPathway

## read genomic matricies
dataFeatures = list(set(mData.retColumns("paradigm/clusterFiles/%s" % (cnvFile))) & set(mData.retColumns("paradigm/clusterFiles/%s" % (expFile))))
dataSamples = list(set(mData.retRows("paradigm/clusterFiles/%s" % (cnvFile))) & set(mData.retRows("paradigm/clusterFiles/%s" % (expFile))))

## load sig_genes.txt
mutData = {}
if os.path.exists("sig_genes.txt"):
    mutData = mData.rCRSData("sig_genes.txt")["p"]

## load include
if os.path.exists("include.samples"):
    dataSamples = list(set(dataSamples) & set(mData.rList("include.samples")))
excludeMap = {}
if os.path.exists("exclude.samples"):
    e = open("exclude.samples", "r")
    for line in e:
        if line.isspace():
            continue
        line.rstrip("\r\n")
        pline = re.split("\t", line)
        excludeMap[pline[0]] = set(re.split(",", pline[1]))

if os.path.exists("include.features"):
    includeFeatures = mData.rList("include.features")
else:
    includeFeatures = None

def writeScripts():
    """creates the R scripts necessary for plotting"""
    signalR = """#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    mutGene = args[1]
    mutFile = args[2]
    wtFile = args[3]

    mutData = read.table(mutFile)
    wtData = read.table(wtFile)
    xMin = min(mutData$V2, wtData$V2)
    if (xMin > 0) {
        xMin = 0
    }
    xMax = max(mutData$V2, wtData$V2)
    if (xMax < 0) {
        xMax = 0
    }
    xMin = 1.4*xMin
    xMax = 1.4*xMax
    
    yMax = 1.1*max(density(mutData$V2)$y, density(wtData$V2)$y)
    
    pdf(paste(mutGene, ".signal.pdf", sep = ""), height = 5, width = 5)
    plot(density(mutData$V2), col = "red", xlim = c(xMin, xMax), ylim = c(0, yMax), main = mutGene, xlab = "")
    par(new = T)
    plot(density(wtData$V2), xlim = c(xMin, xMax), ylim = c(0, yMax), main = "", xlab = "", ylab = "")
    dev.off()
    """
     
    backgroundR = """#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    mutGene = args[1]
    realFile = args[2]
    nullFile = args[3]

    realData = read.table(realFile)$V1
    nullData = read.table(nullFile)$V1
    nullData = nullData[!is.nan(nullData)]
    minVal = min(realData, nullData)
    if (minVal > 0) {
        minVal = 0
    }
    maxVal = max(realData, nullData)
    if (maxVal < 0) {
        maxVal = 0
    }
    minVal = 1.4*minVal
    maxVal = 1.4*maxVal

    pdf(paste(mutGene, ".background.pdf", sep = ""), height = 5, width = 5)
    plot(density(nullData), main = mutGene, xlim = c(minVal, maxVal))
    abline(v = realData[1], col = "red")
    dev.off()
    """
    
    aucR = """#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    mutGene = args[1]
    aucFile = args[2]
    
    aucData = read.delim(aucFile, skip = 1, header = FALSE)
    
    pdf(paste(mutGene, ".auc.pdf", sep = ""), height = 5, width = 5)
    plot(aucData$V1, aucData$V2, pch = "", xlab = "", ylab = "")
    lines(aucData$V1, aucData$V2)
    dev.off()
    """
    
    f = open("signal.R", "w")
    f.write(signalR)
    f.close
    f = open("background.R", "w")
    f.write(backgroundR)
    f.close
    f = open("auc.R", "w")
    f.write(aucR)
    f.close
    system("chmod 755 signal.R background.R auc.R")

class jtData(Target):
    def __init__(self, mutGene, proteinFeatures, directory):
        Target.__init__(self, time=1000)
        self.mutGene = mutGene
        self.proteinFeatures = proteinFeatures
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        random.seed(randSeed)
        
        ## write data
        mData.rwCRSData("full_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), useRows = dataSamples)
        mData.rwCRSData("full_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), useRows = dataSamples)
        mData.rwCRSData("up_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), useRows = dataSamples, useCols = self.proteinFeatures)
        mData.rwCRSData("up_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), useRows = dataSamples, useCols = self.proteinFeatures)
        mData.rwCRSData("down_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), useRows = dataSamples, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))
        mData.rwCRSData("down_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), useRows = dataSamples, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))
        
        ## add perm_ to sample names
        rowMap = {}
        rowFeatures = dataSamples
        rowPermute = ["perm_%s" % (i) for i in rowFeatures]
        for i in range(0, len(rowFeatures)):
            rowMap[rowFeatures[i]] = rowPermute[i]
        
        ## permute features for simulated normals
        rcMap = {}
        for i in rowMap.keys():
            rcMap[rowMap[i]] = {}
            colFeatures = list(set(dataFeatures)-set([self.mutGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[rowMap[i]][colFeatures[j]] = colPermute[j]
        
        ## write simulated normals
        mData.rwCRSData("full_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap)
        mData.rwCRSData("full_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap)
        mData.rwCRSData("up_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("up_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("down_%s" % (cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))
        mData.rwCRSData("down_%s" % (expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))

class jtNData(Target):
    def __init__(self, null, mutGene, proteinFeatures, directory):
        Target.__init__(self, time=1000)
        self.null = null
        self.mutGene = mutGene
        self.proteinFeatures = proteinFeatures
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        random.seed(randSeed+self.null)
        
        ## permute features for nulls
        rcMap = {}
        for i in dataSamples:
            rcMap[i] = {}
            colFeatures = list(set(dataFeatures)-set([self.mutGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[i][colFeatures[j]] = colPermute[j]
        
        ## write data
        mData.rwCRSData("up_N%s_%s" % (self.null, cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("up_N%s_%s" % (self.null, expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("down_N%s_%s" % (self.null, cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))
        mData.rwCRSData("down_N%s_%s" % (self.null, expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))

        ## add perm_ to sample names
        rowMap = {}
        rowFeatures = dataSamples
        rowPermute = ["perm_%s" % (i) for i in rowFeatures]
        for i in range(0, len(rowFeatures)):
            rowMap[rowFeatures[i]] = rowPermute[i]
        
        ## permute features for simulated normals
        rcMap = {}
        for i in rowMap.keys():
            rcMap[rowMap[i]] = {}
            colFeatures = list(set(dataFeatures)-set([self.mutGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[rowMap[i]][colFeatures[j]] = colPermute[j]
                
        ## write simulated normals
        mData.rwCRSData("up_N%s_%s" % (self.null, cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("up_N%s_%s" % (self.null, expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        mData.rwCRSData("down_N%s_%s" % (self.null, cnvFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, cnvFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))
        mData.rwCRSData("down_N%s_%s" % (self.null, expFile), "%s/../paradigm/clusterFiles/%s" % (self.directory, expFile), rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutGene])))

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class branchParams(Target):
    def __init__(self, mutList, sampleMap, directory):
        Target.__init__(self, time=10000)
        self.mutList = mutList
        self.sampleMap = sampleMap
        self.directory = directory
    def run(self):
        for d in distanceParams:
            for t in tstatParams:
                for i in incrParams:
                    for m in methodParams:
                        dirName = "_".join([cohortName, str(d), str(t), str(i), str(m)])
                        system("mkdir %s" % (dirName))
                        self.addChildTarget(branchCV(self.mutList, self.sampleMap, d, t, i, m, "%s/%s" % (self.directory, dirName)))

class branchCV(Target):
    def __init__(self, mutList, sampleMap, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutList = mutList
        self.sampleMap = sampleMap
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        cohortName = re.split("/", self.directory)[-1]
        if os.path.exists("%s/mutpath-%s" % (htmlDir, cohortName)):
            system("rm -rf %s/mutpath-%s" % (htmlDir, cohortName))
        system("mkdir %s/mutpath-%s" % (htmlDir, cohortName))
        if os.path.exists("analysis"):
            system("rm -rf analysis")
        system("mkdir analysis")
        for i in self.mutList:
            system("mkdir analysis/%s" % (i))
            self.addChildTarget(prepareData(i, self.sampleMap[i], list(set(dataSamples) - set(self.sampleMap[i])), self.dParam, self.tParam, self.iParam, self.method, self.directory))
        self.setFollowOnTarget(branchFinal(self.mutList, self.sampleMap, self.dParam, self.tParam, self.iParam, self.method, self.directory))
        
class branchFinal(Target):
    def __init__(self, mutList, sampleMap, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutList = mutList
        self.sampleMap = sampleMap
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        cohortName = re.split("/", self.directory)[-1]
        ## final run
        for i in self.mutList:
            self.addChildTarget(prepareNeighborhood(i, self.sampleMap[i], list(set(dataSamples) - set(self.sampleMap[i])), list(set(dataSamples)), 0, self.dParam, self.tParam, self.iParam, self.method, self.directory))
 
class prepareData(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.nonSamples = nonSamples
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        system("echo Preparing Genomic Data ... >> progress.log")
        
        ## get proteinFeatures from self.dParam and gPathway
        proteinFeatures = []
        for feature in getNeighbors(self.mutGene, self.dParam, gfPathway.interactions):
            if gfPathway.nodes[feature] == "protein":
                proteinFeatures.append(feature)
         
        ## link paradigm requirements
        system("ln -s %s/../paradigm/config.txt config.txt" % (self.directory))
        system("ln -s %s/../paradigm/params.txt params.txt" % (self.directory))
        for i in os.listdir("%s/../paradigm" % (self.directory)):
            if i.endswith(".imap"):
                system("ln -s %s/../paradigm/%s %s" % (self.directory, i, i))
            elif i.endswith(".dogma"):
                system("ln -s %s/../paradigm/%s %s" % (self.directory, i, i))
        
        ## write data
        genData = jtData(self.mutGene, proteinFeatures, self.directory)
        genData.run()
        for i in range(1, nNulls+1):
            genData = jtNData(i, self.mutGene, proteinFeatures, self.directory)
            genData.run()

        ## get paradigm inference
        if useParadigm:
            system("mkdir outputFiles")
            for i in range(len(dataSamples)):
                cmd = "%s -p %s -c config.txt -b ./full_ -o outputFiles/%s_b%s_%s_output.fa -s %s,%s" % (paradigmExec, "%s/../paradigm/clusterFiles/%s" % (self.directory, pathFile), pathFile.rstrip("_pathway.tab"), i, len(dataSamples), i, len(dataSamples))
                self.addChildTarget(jtCmd(cmd, shiftDir))
        self.setFollowOnTarget(prepareMerged(self.mutGene, self.mutSamples, self.nonSamples, self.dParam, self.tParam, self.iParam, self.method, self.directory))
        
class prepareMerged(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.nonSamples = nonSamples
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        
        if useParadigm: 
            system("echo Preparing PARADIGM IPLs ... >> progress.log")
            system("mkdir mergeFiles")
            system("%s outputFiles mergeFiles" % (mswarmExec))
            system("%s bioInt mergeFiles/ merge_merged.pre" % (mmergeExec))
            system("cat merge_merged.pre | sed 's/^[0-9A-Z]*_//' > merge_merged.tab")
            system("cat merge_merged.tab | transpose.pl > merge_merged.transposed.tab")
            system("rm -f merge_merged.pre")
        
        self.setFollowOnTarget(branchFolds(self.mutGene, self.mutSamples, self.nonSamples, self.dParam, self.tParam, self.iParam, self.method, self.directory))

class branchFolds(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.nonSamples = nonSamples
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        system("echo Cross Validation Step ... >> progress.log")
        random.seed(randSeed+123454321)
        
        ## pick samples
        foldSamples = {}
        for r in range(1, rRepeats+1):
            foldSamples[r] = {}
            for f in range(1, mFolds+1):
                foldSamples[r][f] = []
            select_mutSamples = deepcopy(self.mutSamples)
            select_nonSamples = deepcopy(self.nonSamples)
            while len(select_mutSamples)+len(select_nonSamples) > 0:
                for f in range(1, mFolds+1):
                    if len(select_mutSamples) > 0:
                        foldSamples[r][f].append(select_mutSamples.pop(random.randint(0,len(select_mutSamples)-1)))
                    elif len(select_nonSamples) > 0:
                        foldSamples[r][f].append(select_nonSamples.pop(random.randint(0,len(select_nonSamples)-1)))
        
        ## branch folds
        for r in range(1, rRepeats+1):
            for f in range(1, mFolds+1):
                system("mkdir fold%s" % ((r-1)*mFolds+f))
                os.chdir("fold%s" % ((r-1)*mFolds+f))
                mData.wList("include.samples", foldSamples[r][f])
                system("ln -s %s/../paradigm/config.txt config.txt" % (self.directory))
                system("ln -s %s/../paradigm/params.txt params.txt" % (self.directory))
                for i in os.listdir("%s/../paradigm" % (self.directory)):
                    if i.endswith(".imap"):
                        system("ln -s %s/../paradigm/%s %s" % (self.directory, i, i))
                    elif i.endswith(".dogma"):
                        system("ln -s %s/../paradigm/%s %s" % (self.directory, i, i))
                os.chdir("..")
                self.addChildTarget(prepareNeighborhood(self.mutGene, self.mutSamples, self.nonSamples, list(set(dataSamples) - set(foldSamples[r][f])), (r-1)*mFolds+f, self.dParam, self.tParam, self.iParam, self.method, self.directory))

class prepareNeighborhood(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, trainSamples, fold, dParam, tParam, iParam, method, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.nonSamples = nonSamples
        self.trainSamples = trainSamples
        self.fold = fold
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
            os.chdir(shiftDir)
            system("echo Preparing Network ... >> progress.log")
            dataFile = "./up_%s" % (expFile)
            if useParadigm:
                dataFile = "./merge_merged.transposed.tab"
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutGene, self.fold)
            os.chdir(shiftDir)
            dataFile = "../up_%s" % (expFile)
            if useParadigm:
                dataFile = "../merge_merged.transposed.tab"
        
        ## get upstream and downstream
        (uPathway, dPathway, isPass) = selectMutationNeighborhood(self.mutGene, self.mutSamples, dataFile, gPathway, trainSamples = self.trainSamples, tCutoff = self.tParam, tIncrement = self.iParam, maxDistance = self.dParam, method = self.method)
        (uNodes, uInteractions) = (uPathway.nodes, uPathway.interactions)
        (dNodes, dInteractions) = (dPathway.nodes, dPathway.interactions)
        
        ## check thresholds
        if not isPass:
            if self.fold == 0:
                mutGroup = list(set(self.trainSamples) & set(self.mutSamples))
                nonGroup = list(set(self.trainSamples) & set(self.nonSamples))
                f = open("sig.stats", "w")
                f.write("%s\t%s\t%s\t%s\t%s\n" % (self.mutGene, len(nonGroup), len(mutGroup), "NA", "NA"))
                f.close()
                system("cat sig.stats >> ../../sig.stats")
        else:
            self.setFollowOnTarget(runAnalysis(self.mutGene, self.mutSamples, self.nonSamples, self.trainSamples, self.fold, self.directory))

class runAnalysis(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, trainSamples, fold, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.nonSamples = nonSamples
        self.trainSamples = trainSamples
        self.fold = fold
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
            os.chdir(shiftDir)
            system("echo Run PARADIGM ... >> progress.log")
            dataPath = "./"
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutGene, self.fold)
            os.chdir(shiftDir)
            dataPath = "../"
        
        ## run paradigm (observed and nulls)
        self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_" % (dataPath), "%s_upstream.fa" % (self.mutGene)), shiftDir))
        self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_" % (dataPath), "%s_downstream.fa" % (self.mutGene)), shiftDir))
        for i in range(1, nNulls+1):
            self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_N%s_" % (dataPath, i), "N%s_%s_upstream.fa" % (i, self.mutGene)), shiftDir))
            self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_N%s_" % (dataPath, i), "N%s_%s_downstream.fa" % (i, self.mutGene)), shiftDir))
        os.chdir("../..")
        self.setFollowOnTarget(scoreMutation(self.mutGene, self.mutSamples, self.nonSamples, self.trainSamples, self.fold, self.directory))

class scoreMutation(Target):
    def __init__(self, mutGene, mutSamples, nonSamples, trainSamples, fold, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.nonSamples = nonSamples
        self.trainSamples = trainSamples
        self.fold = fold
        self.mutSamples = mutSamples
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
            os.chdir(shiftDir)
            system("echo Evaluate Shifts ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutGene, self.fold)
            os.chdir(shiftDir)
        cohortName = re.split("/", self.directory)[-1]
        
        ## read upstream and downstream pathways
        (uNodes, uInteractions) = rPathway("upstream_pathway.tab")
        (dNodes, dInteractions) = rPathway("downstream_pathway.tab")
        cPathway = combinePathways(Pathway(uNodes, uInteractions), Pathway(dNodes, dInteractions))
        
        ## read paradigm output
        upScore = mData.rPARADIGM("%s_upstream.fa" % (self.mutGene), useRows = uNodes.keys())[1]
        downScore = mData.rPARADIGM("%s_downstream.fa" % (self.mutGene), useRows = dNodes.keys())[1]
        mData.wCRSData("paradigm_up.tab", upScore)
        mData.wCRSData("paradigm_down.tab", downScore)
        n_upScore = {}
        n_downScore = {}
        for i in range(1, nNulls+1):
            n_upScore[i] = mData.rPARADIGM("N%s_%s_upstream.fa" % (i, self.mutGene), useRows = [self.mutGene])[1]
            n_downScore[i] = mData.rPARADIGM("N%s_%s_downstream.fa" % (i, self.mutGene), useRows = [self.mutGene])[1]
        
        ## define sample groups
        mutGroup = self.mutSamples
        nonGroup = self.nonSamples
        permGroup = ["perm_%s" % (i) for i in mutGroup]
        backGroup = []
        for i in mutGroup:
            for j in range(1, nNulls+1):
                backGroup.append("null%s_%s" % (j, i))
        
        ## score shifts      
        shiftScore = {}
        for i in (mutGroup + nonGroup + permGroup):
            shiftScore[i] = downScore[i][self.mutGene] - upScore[i][self.mutGene]
        for i in (mutGroup + nonGroup + permGroup):
            for j in range(1, nNulls+1):
                shiftScore["null%s_%s" % (j, i)] = n_downScore[j][i][self.mutGene] - n_upScore[j][i][self.mutGene]    
        
        ## cross validation
        if self.fold != 0 :
            trainGroup = self.trainSamples
            testGroup = list(set(dataSamples) - set(self.trainSamples))
            mutGroup_tr = list(set(mutGroup) & set(trainGroup))
            mutGroup_te = list(set(mutGroup) & set(testGroup))
            nonGroup_tr = list(set(nonGroup) & set(trainGroup))
            nonGroup_te = list(set(nonGroup) & set(testGroup))
            permGroup_tr = ["perm_%s" % (j) for j in mutGroup_tr]
            permGroup_te = ["perm_%s" % (j) for j in mutGroup_te]
            mData.wList("mut.train", mutGroup_tr)
            mData.wList("non.train", nonGroup_tr)
            mData.wList("perm.train", permGroup_tr)
            (nonMean, nonStd) = mCalculate.mean_std([shiftScore[i] for i in nonGroup_tr])
            (permMean, permStd) = mCalculate.mean_std([shiftScore[i] for i in permGroup_tr])
            
            if usePerm:
                s = open("train.stats", "w")
                for i in mutGroup_tr:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, shiftScore[i]-permMean))
                for i in permGroup_tr:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, shiftScore[i]-permMean))
                s.close()
                s = open("test.stats", "w")
                for i in mutGroup_te:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, shiftScore[i]-permMean))
                for i in permGroup_te:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, shiftScore[i]-permMean))
                s.close()
            
            else:
                s = open("train.stats", "w")
                for i in mutGroup_tr:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, shiftScore[i]-nonMean))
                for i in nonGroup_tr:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, shiftScore[i]-nonMean))
                s.close()
                s = open("test.stats", "w")
                for i in mutGroup_te:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, shiftScore[i]-nonMean))
                for i in nonGroup_te:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, shiftScore[i]-nonMean))
                s.close()
        
        ## final score
        if self.fold == 0:
            ## testing auc
            for r in range(1, rRepeats+1):
                for f in range(1, mFolds+1):
                    system("cat fold%s/test.stats >> test.stats" % ((r-1)*mFolds+f))
                    system("cat fold%s/test.stats | cut -f1,2,3,4 | plotAUC.py -m 2 roc%s.abs.fold%s.txt" % ((r-1)*mFolds+f, self.mutGene, (r-1)*mFolds+f))
            system("cat test.stats | cut -f1,2,3,4 | plotAUC.py -m 2 roc%s.abs.txt" % (self.mutGene))
            system("%s/../auc.R %s roc%s.abs.txt" % (self.directory, self.mutGene, self.mutGene))
            testAUC = mData.retColumns("roc%s.abs.txt" % (self.mutGene))[0]
            
            ## compute statistics for auc and signal plots
            if usePerm:
                s = open("mutation.stats", "w")
                f = open("mut.scores", "w")
                (nonMean, nonStd) = mCalculate.mean_std([shiftScore[i] for i in nonGroup])
                (permMean, permStd) = mCalculate.mean_std([shiftScore[i] for i in permGroup])
                for i in mutGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, permScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-permMean, permScore))
                f.close()
                f = open("non.scores", "w")
                for i in nonGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    # s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, permScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-permMean, permScore))
                f.close()
                f = open("perm.scores", "w")
                for i in permGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, permScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-permMean, permScore))
                f.close()
                s.close()
            else:
                s = open("mutation.stats", "w")
                f = open("mut.scores", "w")
                (nonMean, nonStd) = mCalculate.mean_std([shiftScore[i] for i in nonGroup])
                (permMean, permStd) = mCalculate.mean_std([shiftScore[i] for i in permGroup])
                for i in mutGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 1, nonScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-nonMean, nonScore))
                f.close()
                f = open("non.scores", "w")
                for i in nonGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, nonScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-nonMean, nonScore))
                f.close()
                f = open("perm.scores", "w")
                for i in permGroup:
                    nonScore = (shiftScore[i]-nonMean)/(nonStd+alpha)
                    permScore = (shiftScore[i]-permMean)/(permStd+alpha)
                    # s.write("%s\t%s\t%s\t%s\n" % (self.mutGene, i, 0, nonScore))
                    f.write("%s\t%s\t%s\n" % (i, shiftScore[i]-nonMean, nonScore))
                f.close()
                s.close()
            
            ## compute statistics for background plot
            signalScore = {}
            if usePerm:
                signalScore["real"] = computeSignal(permGroup, mutGroup, shiftScore, method = signalMethod)
                for i in range(1, nNulls+1):
                    signalScore["null%s" % (i)] = computeSignal(["null%s_%s" % (i, j) for j in permGroup], ["null%s_%s" % (i, j) for j in mutGroup], shiftScore, method = signalMethod)

            else:
                signalScore["real"] = computeSignal(nonGroup, mutGroup, shiftScore, method = signalMethod)
                for i in range(1, nNulls+1):
                    signalScore["null%s" % (i)] = computeSignal(["null%s_%s" % (i, j) for j in nonGroup], ["null%s_%s" % (i, j) for j in mutGroup], shiftScore, method = signalMethod)

            f = open("real.scores", "w")
            f.write("%s\n" % (signalScore["real"]))
            f.close()
            f = open("null.scores", "w")
            for i in range(1, nNulls+1):
                f.write("%s\n" % (signalScore["null%s" % (i)]))
            f.close()
            
            ## output signal and background plots
            if usePerm:
                self.addChildTarget(jtCmd("%s/../signal.R %s mut.scores perm.scores" % (self.directory, self.mutGene), shiftDir))
            else:
                self.addChildTarget(jtCmd("%s/../signal.R %s mut.scores non.scores" % (self.directory, self.mutGene), shiftDir))
            self.addChildTarget(jtCmd("%s/../background.R %s real.scores null.scores" % (self.directory, self.mutGene), shiftDir))
            
            ## compute signal strength
            realScores = [signalScore["real"]]
            nullScores = [signalScore[i] for i in list(set(signalScore.keys()) - set(["real"]))]
            (nullMean, nullStd) = mCalculate.mean_std(nullScores)
            strengthScore = (realScores[0]-nullMean)/(nullStd+alpha)
            
            ## write signal strength
            mutVal = "NA"
            if self.mutGene in mutData:
                mutVal = mutData[self.mutGene]
            f = open("sig.stats", "w")
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (self.mutGene, len(nonGroup), len(mutGroup), testAUC, strengthScore, str(mutVal).lstrip("<")))
            f.close()
            system("cat sig.stats >> ../../sig.stats")
            
            ## plot networks
            system("mkdir disc-plot%s" % (self.mutGene))
            system("mkdir disc-plot%s/img" % (self.mutGene))
            proteinFeatures = set()
            gf = open("up.features", "w")
            for i in (set(upScore[upScore.keys()[0]].keys()) - set([self.mutGene])):
                gf.write("%s\n" % (i))
                proteinFeatures.update([i])
            gf.close()
            gf = open("down.features", "w")
            for i in (set(downScore[downScore.keys()[0]].keys()) - set([self.mutGene])):
                gf.write("%s\n" % (i))
                proteinFeatures.update([i])
            gf.close()
            gf = open("mut.features", "w")
            gf.write("%s\n" % (self.mutGene))
            gf.close()
            if usePerm:
                sf = open("include.samples", "w")
                for i in (mutGroup + nonGroup + permGroup):
                    sf.write("%s\n" % (i))
                sf.close()
            else:
                sf = open("include.samples", "w")
                for i in (mutGroup + nonGroup):
                    sf.write("%s\n" % (i))
                sf.close()
            mf = open("mut.circle", "w")
            mf.write("id\t%s\n" % ("\t".join(mutGroup + nonGroup + permGroup)))
            mf.write("*")
            for i in (mutGroup + nonGroup + permGroup):
                if i in mutGroup:
                    mf.write("\t1")
                elif i in nonGroup:
                    mf.write("\t0")
                else:
                    mf.write("\t0.5")
            mf.write("\n")
            mf.close()
            mf = open("shift.circle", "w")
            mf.write("id\t%s\n" % ("\t".join(mutGroup + nonGroup + permGroup)))
            mf.write("*")
            for i in (mutGroup + nonGroup + permGroup):
                mf.write("\t%s" % (shiftScore[i]))
            mf.write("\n")
            mf.close()
            system("cat up_%s | transpose.pl > exp_t.tab" % (expFile))
            system("%s -s \"-0.5\" exp_t.tab exp.tab" % (medianExec))
            system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f mut.features disc-plot%s/img/ mut.circle exp.tab paradigm_up.tab paradigm_down.tab shift.circle" % (circleExec, self.mutGene, self.mutGene))
            system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f up.features disc-plot%s/img/ mut.circle exp.tab paradigm_up.tab" % (circleExec, self.mutGene, self.mutGene))
            system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f down.features disc-plot%s/img/ mut.circle exp.tab paradigm_down.tab" % (circleExec, self.mutGene, self.mutGene))
            scoreMap = {}
            scoreMap["disc-plot%s" % (self.mutGene)] = {}
            for node in gPathway.nodes.keys():
                if node == self.mutGene:
                    scoreMap["disc-plot%s" % (self.mutGene)][node] = 10
                else:
                    scoreMap["disc-plot%s" % (self.mutGene)][node] = 7
            wSIF("disc-plot%s/disc-plot%s.sif" % (self.mutGene, self.mutGene), cPathway.interactions)
            #wSIF("disc-plot%s/disc-plot%s.sif_with_complexes" % (self.mutGene, self.mutGene), cPathway.interactions)
            #cfPathway = flattenPathway(cPathway)
            #wSIF("disc-plot%s/disc-plot%s.sif" % (self.mutGene, self.mutGene), cfPathway.interactions)
            wNodeAttributes(gPathway.nodes, scoreMap = scoreMap, directory = "./")
            
            ## cleanup
            self.setFollowOnTarget(cleanUp(self.mutGene, self.mutSamples, self.directory))

class cleanUp(Target):
    def __init__(self, mutGene, mutSamples, directory):
        Target.__init__(self, time=10000)
        self.mutGene = mutGene
        self.mutSamples = mutSamples
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutGene)
        os.chdir(shiftDir)
        system("echo Plotting CytoscapeWeb ... >> progress.log")
        cohortName = re.split("/", self.directory)[-1]
        
        ## cytoscape-web
        system("%s disc-plot%s %s/mutpath-%s" % (cytowebExec, self.mutGene, htmlDir, cohortName))
        system("cat mutation.stats | cut -f2- | cap.pl \"id,MutStatus,Score\" | tab2html.py %s/mutpath-%s/table-%s.html" % (htmlDir, cohortName, self.mutGene))
        
        ## cleanup
        system("rm -f include.samples up.features down.features mut.circle exp.tab exp_t.tab nullexp.tab nullexp_t.tab paradigm_down.tab paradigm_up.tab nullparadigm_down.tab nullparadigm_up.tab shift.circle")
        system("rm -f up_* down_*")
        system("rm -f *.fa")
        
        ## update
        os.chdir(htmlDir)
        system("./update-html.py")
        
def main():
    ## parse arguments
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help="Add as a child of jobFile rather " +
                      "than making a new jobTree")
    options, args = parser.parse_args()
    print "Using Batch System '" + options.batchSystem + "'"
    ## clean
    if len(args) == 1:
        if args[0] == "clean":
            print "rm -rf %s* %s/mutpath-%s* .jobTree" % (cohortName, htmlDir, cohortName)
            os.system("rm -rf %s* %s/mutpath-%s* .jobTree"% (cohortName, htmlDir, cohortName))
            sys.exit(0)
    
    ## run
    assert len(args) == 0
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    mutGenes = []
    mutSamples = {}
    f = open("mutationSummary.tab", "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split("\t", line)
        gene = pline[0]
        samples = list(set(re.split(",", pline[2])) & set(dataSamples))
        if gene in gPathway.nodes:
            if len(samples) >= mutThresh:
                if includeFeatures != None:
                    if gene not in includeFeatures:
                        continue
                mutGenes.append(gene)
                mutSamples[gene] = deepcopy(samples)
    f.close()
    s = Stack(branchParams(mutGenes, mutSamples, os.getcwd()))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = "./.jobTree"
        
        failed = s.startJobTree(options)
        if failed:
            print ("%d jobs failed" % failed)
        else:
            logger.info("Run complete!")
            system("rm -rf .lastjobTree")
            system("mv .jobTree .lastjobTree")

if __name__ == "__main__":
    from paradigmSHIFT import *
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory already exists"
    main()
