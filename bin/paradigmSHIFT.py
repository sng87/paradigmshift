#!/usr/bin/env python
"""
paradigmSHIFT.py
"""
## Written By: Sam Ng
import math, os, sys, random, re
from copy import deepcopy

from mSHIFT import *

from optparse import OptionParser
from jobTree.src.bioio import logger
from jobTree.src.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## randomize seed
if not os.path.exists("seed.log"):
    randSeed = random.randint(0,9999999999999999)
    f = open("seed.log", "w")
    f.write("%s\n" % (randSeed))
    f.close()
else:
    f = open("seed.log", "r")
    randSeed = int(f.readline().rstrip("\r\n"))
    f.close()

## executables and directories
paradigmExec = "paradigm"
circleExec = "circlePlot.py"
htmlDir = "/inside/grotto/users/sng/.html"

## default variables
paradigmPublic = False            ## Use public binary (does not support sample splitting)
useGreedy = False                 ## Use greedy network search (not implemented)
useFlattened = True               ## Flattens complexes to prevent long complex chains
maxFeatures = 10                  ## Maximum number of features to in PARADIGM Shift run
mutationThreshold = 5             ## Minimum mutations in cohort needed to consider a gene
nNulls = 30                       ## Number of nulls to run for background model

def writeScripts():
    """creates the R scripts necessary for plotting"""
    msepR = """#!/usr/bin/env Rscript
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
    
    pdf(paste(mutGene, ".msep.pdf", sep = ""), height = 5, width = 5)
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
    
    f = open("msep.R", "w")
    f.write(msepR)
    f.close
    f = open("background.R", "w")
    f.write(backgroundR)
    f.close
    system("chmod 755 msep.R background.R")

class jtData(Target):
    def __init__(self, mutatedGene, dataSamples, proteinFeatures, dataFeatures, dataMap, directory):
        Target.__init__(self, time=1000)
        self.mutatedGene = mutatedGene
        self.dataSamples = dataSamples
        self.proteinFeatures = proteinFeatures
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        random.seed(randSeed)
        
        ## get files
        cnvFile = self.dataMap["cnv"]
        expFile = self.dataMap["exp"]
        
        ## write data
        rwCRSData("full_%s" % (re.split("/", cnvFile)[-1]), cnvFile, useRows = self.dataSamples)
        rwCRSData("full_%s" % (re.split("/", expFile)[-1]), expFile, useRows = self.dataSamples)
        rwCRSData("up_%s" % (re.split("/", cnvFile)[-1]), cnvFile, useRows = self.dataSamples, useCols = self.proteinFeatures)
        rwCRSData("up_%s" % (re.split("/", expFile)[-1]), expFile, useRows = self.dataSamples, useCols = self.proteinFeatures)
        rwCRSData("down_%s" % (re.split("/", cnvFile)[-1]), cnvFile, useRows = self.dataSamples, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))
        rwCRSData("down_%s" % (re.split("/", expFile)[-1]), expFile, useRows = self.dataSamples, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))
        
        ## add perm_ to sample names
        rowMap = {}
        rowFeatures = self.dataSamples
        rowPermute = ["perm_%s" % (i) for i in rowFeatures]
        for i in range(0, len(rowFeatures)):
            rowMap[rowFeatures[i]] = rowPermute[i]
        
        ## permute features for simulated normals
        rcMap = {}
        for i in rowMap.keys():
            rcMap[rowMap[i]] = {}
            colFeatures = list(set(self.dataFeatures)-set([self.mutatedGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[rowMap[i]][colFeatures[j]] = colPermute[j]
        
        ## write simulated normals
        rwCRSData("full_%s" % (re.split("/", cnvFile)[-1]), cnvFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap)
        rwCRSData("full_%s" % (re.split("/", expFile)[-1]), expFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap)
        rwCRSData("up_%s" % (re.split("/", cnvFile)[-1]), cnvFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("up_%s" % (re.split("/", expFile)[-1]), expFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("down_%s" % (re.split("/", cnvFile)[-1]), cnvFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))
        rwCRSData("down_%s" % (re.split("/", expFile)[-1]), expFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))

class jtNData(Target):
    def __init__(self, null, mutatedGene, dataSamples, proteinFeatures, dataFeatures, dataMap, directory):
        Target.__init__(self, time=1000)
        self.null = null
        self.mutatedGene = mutatedGene
        self.dataSamples = dataSamples
        self.proteinFeatures = proteinFeatures
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        random.seed(randSeed+self.null)
        
        ## get files
        cnvFile = self.dataMap["cnv"]
        expFile = self.dataMap["exp"]
        
        ## permute features for nulls
        rcMap = {}
        for i in self.dataSamples:
            rcMap[i] = {}
            colFeatures = list(set(self.dataFeatures)-set([self.mutatedGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[i][colFeatures[j]] = colPermute[j]
        
        ## write data
        rwCRSData("up_N%s_%s" % (self.null, re.split("/", cnvFile)[-1]), cnvFile, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("up_N%s_%s" % (self.null, re.split("/", expFile)[-1]), expFile, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("down_N%s_%s" % (self.null, re.split("/", cnvFile)[-1]), cnvFile, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))
        rwCRSData("down_N%s_%s" % (self.null, re.split("/", expFile)[-1]), expFile, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))

        ## add perm_ to sample names
        rowMap = {}
        rowFeatures = self.dataSamples
        rowPermute = ["perm_%s" % (i) for i in rowFeatures]
        for i in range(0, len(rowFeatures)):
            rowMap[rowFeatures[i]] = rowPermute[i]
        
        ## permute features for simulated normals
        rcMap = {}
        for i in rowMap.keys():
            rcMap[rowMap[i]] = {}
            colFeatures = list(set(self.dataFeatures)-set([self.mutatedGene]))
            colPermute = random.sample(colFeatures, len(colFeatures))
            for j in range(0, len(colFeatures)):
                rcMap[rowMap[i]][colFeatures[j]] = colPermute[j]
                
        ## write simulated normals
        rwCRSData("up_N%s_%s" % (self.null, re.split("/", cnvFile)[-1]), cnvFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("up_N%s_%s" % (self.null, re.split("/", expFile)[-1]), expFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = self.proteinFeatures)
        rwCRSData("down_N%s_%s" % (self.null, re.split("/", cnvFile)[-1]), cnvFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))
        rwCRSData("down_N%s_%s" % (self.null, re.split("/", expFile)[-1]), expFile, rowMap = rowMap, useRows = rcMap.keys(), rcMap = rcMap, useCols = list(set(self.proteinFeatures)-set([self.mutatedGene])))

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class branchGenes(Target):
    def __init__(self, dataSamples, dataFeatures, dataMap, mutationMap, gPathway, paradigmDir, 
                 paramMap, foldMap, directory):
        Target.__init__(self, time=10000)
        self.dataSamples = dataSamples
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.mutationMap = mutationMap
        self.gPathway = gPathway
        self.paradigmDir = paradigmDir
        self.paramMap = paramMap
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## branch genes
        htmlFeatures = []
        if not os.path.exists("analysis"):
            system("mkdir analysis")
        for mutatedGene in self.mutationMap.keys():
            if not os.path.exists("analysis/%s" % (mutatedGene)):
                system("mkdir analysis/%s" % (mutatedGene))
                htmlFeatures.append(mutatedGene)
                self.addChildTarget(branchFolds(mutatedGene, self.mutationMap[mutatedGene], 
                                            self.dataSamples, self.dataFeatures, self.dataMap, 
                                            self.gPathway, self.paradigmDir, self.paramMap, 
                                            self.foldMap, self.directory))
        if os.path.exists(htmlDir):
            self.setFollowOnTarget(pshiftReport(htmlFeatures, "%s/%s" % (htmlDir, self.paramMap["cohortName"]), self.directory))

class branchFolds(Target):
    def __init__(self, mutatedGene, mutatedSamples, dataSamples, dataFeatures, dataMap, gPathway, 
                 paradigmDir, paramMap, foldMap, directory):
        Target.__init__(self, time=10000)
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.gPathway = gPathway
        self.paradigmDir = paradigmDir
        self.paramMap = paramMap
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        random.seed(randSeed+123454321)
        
        ## prepare htmlDir
        
        ## get proteinFeatures based on maxDist  and gPathway
        system("echo Preparing Genomic Data ... >> progress.log")
        proteinFeatures = []
        for feature in getNeighbors(self.mutatedGene, max(self.paramMap["dist"]) + 1, self.gPathway.interactions):
            if self.gPathway.nodes[feature] == "protein":
                proteinFeatures.append(feature)
        
        ## write data
        genData = jtData(self.mutatedGene, self.dataSamples, proteinFeatures, self.dataFeatures, 
                         self.dataMap, self.directory)
        genData.run()
        for null in range(1, nNulls+1):
            genData = jtNData(null, self.mutatedGene, self.dataSamples, proteinFeatures, 
                              self.dataFeatures, self.dataMap, self.directory)
            genData.run()
        
        ## pick samples
        system("echo Branching Folds and Params... >> progress.log")
        rRepeats = len(self.foldMap.keys())
        mFolds = len(self.foldMap[1].keys())
        if self.foldMap[1][1] is None:
            mutSamples = self.mutatedSamples
            nonSamples = list(set(self.dataSamples) - set(self.mutatedSamples))
            foldSamples = {}
            for r in range(1, rRepeats+1):
                foldSamples[r] = {}
                for f in range(1, mFolds+1):
                    foldSamples[r][f] = []
                select_mutSamples = deepcopy(mutSamples)
                select_nonSamples = deepcopy(nonSamples)
                while len(select_mutSamples)+len(select_nonSamples) > 0:
                    for f in range(1, mFolds+1):
                        if len(select_mutSamples) > 0:
                            foldSamples[r][f].append(select_mutSamples.pop(random.randint(0, 
                                                                           len(select_mutSamples)-1)))
                        elif len(select_nonSamples) > 0:
                            foldSamples[r][f].append(select_nonSamples.pop(random.randint(0, 
                                                                           len(select_nonSamples)-1)))
        else:
            foldSamples = deepcopy(self.foldMap)
        
        ## branch folds
        for r in range(1, rRepeats+1):
            for f in range(1, mFolds+1):
                fold = (r-1)*mFolds+f
                system("mkdir fold%s" % (fold))
                self.addChildTarget(branchParams(fold, self.mutatedGene, self.mutatedSamples, 
                                                 self.dataSamples, foldSamples[r][f], 
                                                 self.dataFeatures, self.dataMap, self.gPathway,
                                                 self.paradigmDir, self.paramMap, self.foldMap, 
                                                 self.directory))
        
        ## run final
        self.setFollowOnTarget(branchParams(0, self.mutatedGene, self.mutatedSamples, 
                                            self.dataSamples, self.dataSamples, self.dataFeatures, 
                                            self.dataMap, self.gPathway, self.paradigmDir, 
                                            self.paramMap, self.foldMap, self.directory))

class branchParams(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, dataFeatures, 
                 dataMap, gPathway, paradigmDir, paramMap, foldMap, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.gPathway = gPathway
        self.paradigmDir = paradigmDir
        self.paramMap = paramMap
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
            os.chdir(shiftDir)
            system("echo Building Final Model ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutatedGene, self.fold)
            os.chdir(shiftDir)
        
        ## average cross validation aucs for final run
        rRepeats = len(self.foldMap.keys())
        mFolds = len(self.foldMap[1].keys())
        if self.fold == 0:
            aucList = []
            aucLines = []
            for r in range(1, rRepeats+1):
                for f in range(1, mFolds+1):
                    fold = (r-1)*mFolds+f
                    f = open("fold%s/auc.stat" % (fold), "r")
                    line = f.readline()
                    f.close()
                    (auc_tr, auc_te, mparams) = re.split("\t", line.rstrip("\r\n"))
                    aucList.append(auc_te)
                    aucLines.append("%s\t%s\t%s\t%s\n" % (fold, auc_tr, auc_te, mparams))
            o = open("avgAUC.tab", "w")
            o.write("> %s\tavgAUC:%s\n" % (self.mutatedGene, mean(aucList)))
            o.write("# fold\ttrain\ttest\tparams\n")
            for line in aucLines:
                o.write(line)
            o.close()
        
        ## branch params
        for dist in self.paramMap["dist"]:
            for thresh in self.paramMap["thresh"]:
                for inc in self.paramMap["inc"]:
                    for method in self.paramMap["method"]:
                        system("mkdir param_%s_%s_%s_%s" % (dist, thresh, inc, method))
                        os.chdir("param_%s_%s_%s_%s" % (dist, thresh, inc, method))
                        system("cp %s/config.txt config.txt" % (self.paradigmDir))
                        system("cp %s/params.txt params.txt" % (self.paradigmDir))
                        for file in os.listdir(self.paradigmDir):
                            if file.endswith(".imap"):
                                system("cp %s/%s %s" % (self.paradigmDir, file, file))
                            elif file.endswith(".dogma"):
                                system("cp %s/%s %s" % (self.paradigmDir, file, file))
                        os.chdir("..")
                        self.addChildTarget(prepareNeighborhood(self.fold, self.mutatedGene, 
                                                    self.mutatedSamples, self.dataSamples, 
                                                    self.trainSamples, dist, thresh, inc, method,
                                                    self.dataMap, self.gPathway, 
                                                    self.directory))
        
        ## evaluate models
        self.setFollowOnTarget(evaluateParams(self.fold, self.mutatedGene, self.mutatedSamples, 
                                              self.dataSamples, self.trainSamples, 
                                              self.dataFeatures, self.dataMap, self.gPathway, 
                                              self.paramMap, self.directory))
        
class prepareNeighborhood(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, dParam,
                 tParam, iParam, method, dataMap, gPathway, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.dataMap = dataMap
        self.gPathway = gPathway
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s/param_%s_%s_%s_%s" % (self.directory, self.mutatedGene, 
                                                             self.dParam, self.tParam, self.iParam, 
                                                             self.method)
            os.chdir(shiftDir)
            system("echo Preparing Network ... >> progress.log")
            dataFile = "../up_%s" % (re.split("/", self.dataMap["exp"])[-1])
        else:
            shiftDir = "%s/analysis/%s/fold%s/param_%s_%s_%s_%s" % (self.directory,
                                                                    self.mutatedGene, self.fold, 
                                                                    self.dParam, self.tParam, 
                                                                    self.iParam, self.method)
            os.chdir(shiftDir)
            dataFile = "../../up_%s" % (re.split("/", self.dataMap["exp"])[-1])
        
        ## get upstream and downstream
        system("echo Selecting Mutational Network ... >> progress.log")
        (uPathway, dPathway, isPass) = selectMutationNeighborhood(self.mutatedGene, 
                                                                  self.mutatedSamples, dataFile, 
                                                                  self.gPathway, 
                                                                  trainSamples = self.trainSamples, 
                                                                  tCutoff = self.tParam, 
                                                                  tIncrement = self.iParam, 
                                                                  maxDistance = self.dParam, 
                                                                  method = self.method)
        wPathway("upstream_pathway.tab", uPathway.nodes, uPathway.interactions)
        wPathway("downstream_pathway.tab", dPathway.nodes, dPathway.interactions)
        
        ## check thresholds
        if not isPass:
            system("echo Stop ... >> progress.log")
            f = open("auc.stat", "w")
            f.write("NA\tNA\n")
            f.close()
        else:
            if useGreedy:
                self.setFollowOnTarget(runGreedy(self.fold, self.mutatedGene, self.mutatedSamples, 
                                                 self.dataSamples, self.trainSamples, uPathway, 
                                                 dPathway, self.directory, shiftDir))
            else:
                if self.fold == 0:
                    dataPath = "../"
                else:
                    dataPath = "../../"
                self.setFollowOnTarget(runPARADIGM(self.fold, self.mutatedGene, self.mutatedSamples, 
                                                self.dataSamples, self.trainSamples, uPathway, 
                                                dPathway, dataPath, shiftDir))

class runGreedy(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, 
                 dPathway, root, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.root = root
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## assert links.tab exists
        assert os.path.exists("%s/links.tab" % (self.root))
        (lNodes, lInteractions) = rPathway("%s/links.tab" % (self.root))
        lPathway = Pathway(lNodes, lInteractions)
        
        ## determine possible additions
        upstreamAdditions = []
        downstreamAdditions = []
        for source in lInteractions.keys():
            for target in lInteractions[source].keys():
                if target in self.uPathway.nodes:
                    if lInteractions[source][target].startswith("-a"):
                        if source in self.uPathway.interactions:
                            if target not in self.uPathway.interactions[source]:
                                upstreamAdditions.append( (source, target, lInteractions[source][target]) )
                        else:
                            upstreamAdditions.append( (source, target, lInteractions[source][target]) )
                ### update this to include downstream chains
                if source in self.dPathway.nodes:
                    if lInteractions[source][target].startswith("-t"):
                        if source in self.dPathway.interactions:
                            if target not in self.uPathway.interactions[source]:
                                downstreamAdditions.append( (source, target, lInteractions[source][target]) )
                        else:
                            downstreamAdditions.append( (source, target, lInteractions[source][target]) )
        additionMap = {}
        for index in range(upstreamAdditions):
            additionMap["u%s" % (index)] = upstreamAdditions[index]
        for index in range(downstreamAdditions):
            additionMap["d%s" % (index)] = downstreamAdditions[index]
        
        ## get baseline PARADIGM-SHIFT auc
        system("mkdir iter_0")
        wPathway("iter_0/upstream_pathway.tab", self.uPathway.nodes, self.uPathway.interactions)
        wPathway("iter_0/downstream_pathway.tab", self.dPathway.nodes, self.dPathway.interactions)
        system("cp config.txt iter_0/config.txt")
        system("cp params.txt iter_0/params.txt")
        for file in os.listdir("."):
            if file.endswith(".imap"):
                system("cp %s iter_0/%s" % (file, file))
            elif file.endswith(".dogma"):
                system("cp %s iter_0/%s" % (file, file))
        if self.fold == 0:
            dataPath = "../../"
        else:
            dataPath = "../../../"
        self.addChildTarget(runPARADIGM(self.fold, self.mutatedGene, self.mutatedSamples, 
                                        self.dataSamples, self.trainSamples, self.uPathway, 
                                        self.dPathway, dataPath, "%s/iter_0" % (self.directory)))
        self.setFollowOnTarget(greedyBranchSearch(self.fold, self.mutatedGene, self.mutatedSamples, 
                                        self.dataSamples, self.trainSamples, self.uPathway, 
                                        self.dPathway, 1, 5, [], additionMap, lPathway, self.directory))
        
class greedyBranchSearch(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, 
                 dPathway, iter, maxiter, additions, additionMap, links, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.iter = iter
        self.maxiter = maxiter
        self.additions = additions
        self.additionMap = additionMap
        self.links = links
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## get current base
        baseUp = deepcopy(self.uPathway)
        baseDown = deepcopy(self.dPathway)
        for addition in self.additions:
            ### update this to include chains
            if addition.startswith("u"):
                if self.additionMap[addition][0] not in baseUp.nodes:
                    baseUp.nodes[self.additionMap[addition][0]] = self.links.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][1] not in baseUp.nodes:
                    baseUp.nodes[self.additionMap[addition][1]] = self.links.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in baseUp.interactions:
                    baseUp.interactions[self.additionMap[addition][0]] = {}
                baseUp.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            else:
                if self.additionMap[addition][0] not in baseDown.nodes:
                    baseDown.nodes[self.additionMap[addition][0]] = self.links.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][1] not in baseDown.nodes:
                    baseDown.nodes[self.additionMap[addition][1]] = self.links.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in baseDown.interactions:
                    baseDown.interactions[self.additionMap[addition][0]] = {}
                baseDown.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
        
        ## for each addition
        for addition in self.additionMap.keys():
            if addition in self.additions:
                continue
            currentUp = deepcopy(baseUp)
            currentDown = deepcopy(baseDown)
            ### update this to include chains
            if addition.startswith("u"):
                if self.additionMap[addition][0] not in currentUp.nodes:
                    currentUp.nodes[self.additionMap[addition][0]] = self.links.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][1] not in currentUp.nodes:
                    currentUp.nodes[self.additionMap[addition][1]] = self.links.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in currentUp.interactions:
                    currentUp.interactions[self.additionMap[addition][0]] = {}
                currentUp.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            else:
                if self.additionMap[addition][0] not in currentDown.nodes:
                    currentDown.nodes[self.additionMap[addition][0]] = self.links.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][1] not in currentDown.nodes:
                    currentDown.nodes[self.additionMap[addition][1]] = self.links.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in currentDown.interactions:
                    currentDown.interactions[self.additionMap[addition][0]] = {}
                currentDown.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            system("mkdir iter_%s-%s" % (self.iter, addition))
            wPathway("iter_%s-%s/upstream_pathway.tab" % (self.iter, addition), currentUp.nodes, currentUp.interactions)
            wPathway("iter_%s-%s/downstream_pathway.tab" % (self.iter, addition), currentDown.nodes, currentDown.interactions)
            system("cp config.txt iter_%s-%s/config.txt" % (self.iter, addition))
            system("cp params.txt iter_%s-%s/params.txt" % (self.iter, addition))
            for file in os.listdir("."):
                if file.endswith(".imap"):
                    system("cp %s iter_%s-%s/%s" % (file, self.iter, addition, file))
                elif file.endswith(".dogma"):
                    system("cp %s iter_%s-%s/%s" % (file, self.iter, addition, file))
            if self.fold == 0:
                dataPath = "../../"
            else:
                dataPath = "../../../"
            self.addChildTarget(runPARADIGM(self.fold, self.mutatedGene, self.mutatedSamples, 
                                            self.dataSamples, self.trainSamples, self.uPathway, 
                                            self.dPathway, dataPath, 
                                            "%s/iter_%s-%s" % (self.directory, self.iter, addition)))
        self.setFollowOnTarget(greedySelectOptimum(self.fold, self.mutatedGene, self.mutatedSamples, 
                                                   self.dataSamples, self.trainSamples, 
                                                   self.uPathway, self.dPathway, self.iter, 
                                                   self.maxiter, self.additions, self.additionMap,
                                                   self.links, self.directory))      
            
class greedySelectOptimum(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, 
                 dPathway, iter, maxiter, additions, additionMap, links, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.iter = iter
        self.maxiter = maxiter
        self.additions = additions
        self.additionMap = additionMap
        self.links = links
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        aucDiffs = {}
        f = open("iter_%s/auc.stat" % (self.iter-1), "r")
        lastAUC = float(f.readline().rstrip("\r\n").split("\t")[0])
        f.close()
        for addition in self.additionMap.keys():
            if os.path.exists("iter_%s-%s/auc.stat" % (self.iter, addition)):
                f = open("iter_%s-%s/auc.stat" % (self.iter, addition), "r")
                currentAUC = float(f.readline().rstrip("\r\n").split("\t")[0])
                f.close()
                aucDiffs[addition] = currentAUC - lastAUC
        rankedAdditions = aucDiffs.keys()
        rankedAdditions.sort(lambda x, y: cmp(aucDiffs[y], aucDiffs[x]))
        if aucDiffs[rankedAdditions[0]] > 0.0:
            updatedAdditions = deepcopy(self.additions)
            updatedAdditions.append(rankedAdditions[0])
            system("ln -s iter_%s-%s iter_%s" % (self.iter, rankedAdditions[0], self.iter))
            if self.iter < self.maxiter:
                self.setFollowOnTarget(greedyBranchSearch(self.fold, self.mutatedGene, 
                                                          self.mutatedSamples, self.dataSamples, 
                                                          self.trainSamples, self.uPathway, 
                                                          self.dPathway, self.iter+1, self.maxiter, 
                                                          updatedAdditions, self.additionMap, 
                                                          self.links, self.directory))
            else:
                for file in ["shift.tab", "include.samples", "up.features", "down.features", 
                             "mut.features", "mut.circle", "shift.circle", "pshift.train.tab",
                             "mut.train.scores", "non.train.scores", "pshift.test.tab",
                             "mut.test.scores", "non.test.scores", "auc.stat", "real.scores",
                             "null.scores", "sig.tab"]:
                    if os.path.exists("iter_%s/%s" % (self.iter, file)):
                        system("cp iter_%s/%s %s" % (self.iter, file, file))        
        else:
            for file in ["shift.tab", "include.samples", "up.features", "down.features", 
                         "mut.features", "mut.circle", "shift.circle", "pshift.train.tab",
                         "mut.train.scores", "non.train.scores", "pshift.test.tab",
                         "mut.test.scores", "non.test.scores", "auc.stat", "real.scores",
                         "null.scores", "sig.tab"]:
                if os.path.exists("iter_%s/%s" % (self.iter-1, file)):
                    system("cp iter_%s/%s %s" % (self.iter-1, file, file))

class runPARADIGM(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, 
                 dPathway, dataPath, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.dataPath = dataPath
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## run paradigm (observed and nulls)
        system("echo Running PARADIGM inference ... >> progress.log")
        if not paradigmPublic:
			self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_" % (self.dataPath), "%s_upstream.fa" % (self.mutatedGene)), self.directory))
			self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_" % (self.dataPath), "%s_downstream.fa" % (self.mutatedGene)), self.directory))
			for null in range(1, nNulls+1):
				self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_N%s_" % (self.dataPath, null), "N%s_%s_upstream.fa" % (null, self.mutatedGene)), self.directory))
				self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_N%s_" % (self.dataPath, null), "N%s_%s_downstream.fa" % (null, self.mutatedGene)), self.directory))
        else:
            system("mkdir outputFiles")
            for b in range(len(self.dataSamples)):
                self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/up_" % (self.dataPath), "outputFiles/%s_upstream_b%s_%s.fa" % (self.mutatedGene, b, len(self.dataSamples)), b, len(self.dataSamples)), self.directory))
                self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/down_" % (self.dataPath), "outputFiles/%s_downstream_b%s_%s.fa" % (self.mutatedGene, b, len(self.dataSamples)), b, len(self.dataSamples)), self.directory))
                for null in range(1, nNulls+1):
                    self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/up_N%s_" % (self.dataPath, null), "outputFiles/N%s_%s_upstream_b%s_%s.fa" % (null, self.mutatedGene, b, len(self.dataSamples)), b, len(self.dataSamples)), self.directory))
                    self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/down_N%s_" % (self.dataPath, null), "outputFiles/N%s_%s_downstream_b%s_%s.fa" % (null, self.mutatedGene, b, len(self.dataSamples)), b, len(self.dataSamples)), self.directory))
        self.setFollowOnTarget(evaluateCV(self.fold, self.mutatedGene, self.mutatedSamples, 
                                          self.dataSamples, self.trainSamples, self.uPathway, 
                                          self.dPathway, self.directory))

class evaluateCV(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, 
                 dPathway, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        if paradigmPublic:
            for b in range(len(self.dataSamples)):
                system("cat outputFiles/%s_upstream_b%s_%s.fa >> %s_upstream.fa" % (self.mutatedGene, b, len(self.dataSamples), self.mutatedGene))
                system("cat outputFiles/%s_downstream_b%s_%s.fa >> %s_downstream.fa" % (self.mutatedGene, b, len(self.dataSamples), self.mutatedGene))
                for null in range(1, nNulls+1):
                    system("cat outputFiles/N%s_%s_upstream_b%s_%s.fa >> N%s_%s_upstream.fa" % (null, self.mutatedGene, b, len(self.dataSamples), null, self.mutatedGene))
                    system("cat outputFiles/N%s_%s_downstream_b%s_%s.fa >> N%s_%s_downstream.fa" % (null, self.mutatedGene, b, len(self.dataSamples), null, self.mutatedGene))
            system("rm -rf outputFiles")
        
        shiftCV(self.mutatedGene, self.mutatedSamples, self.dataSamples, self.trainSamples, 
                self.uPathway, self.dPathway, nNulls = nNulls)

class evaluateParams(Target):
    def __init__(self, fold, mutatedGene, mutatedSamples, dataSamples, trainSamples, dataFeatures, 
                 dataMap, gPathway, paramMap, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.trainSamples = trainSamples
        self.dataFeatures = dataFeatures
        self.dataMap = dataMap
        self.gPathway = gPathway
        self.paramMap = paramMap
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
            os.chdir(shiftDir)
            system("echo Identifying Best Param by Training Validation ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutatedGene, self.fold)
            os.chdir(shiftDir)
        
        ## report AUCs for model with best training validation
        topAUC_tr = 0
        topAUC_te = 0
        topAUC_params = None
        for dist in self.paramMap["dist"]:
            for thresh in self.paramMap["thresh"]:
                for inc in self.paramMap["inc"]:
                    for method in self.paramMap["method"]:
                        f = open("param_%s_%s_%s_%s/auc.stat" % (dist, thresh, inc, method), "r")
                        line = f.readline()
                        f.close()
                        (auc_tr, auc_te) = re.split("\t", line.rstrip("\r\n"))
                        try:
                            auc_tr = float(auc_tr)
                        except ValueError:
                            continue
                        try:
                            auc_te = float(auc_te)
                        except ValueError:
                            pass
                        if auc_tr > topAUC_tr:
                            topAUC_tr = auc_tr
                            topAUC_te = auc_te
                            topAUC_params = [str(dist), str(thresh), str(inc), str(method)]
        if self.fold != 0:
            f = open("auc.stat", "w")
            if topAUC_tr == 0:
                f.write("NA\tNA\tNA\n")
            else:
                f.write("%s\t%s\t%s\n" % (topAUC_tr, topAUC_te, ",".join(topAUC_params)))
            f.close()
        else:
            if topAUC_params is not None:
                self.setFollowOnTarget(mutationSHIFT(topAUC_params, self.mutatedGene, 
                                                     self.mutatedSamples, self.dataSamples, 
                                                     self.dataMap, self.gPathway, self.directory))
            
class mutationSHIFT(Target):
    def __init__(self, topAUC_params, mutatedGene, mutatedSamples, dataSamples, 
                 dataMap, gPathway, directory):
        Target.__init__(self, time=10000)
        self.topAUC_params = topAUC_params
        self.mutatedGene = mutatedGene
        self.mutatedSamples = mutatedSamples
        self.dataSamples = dataSamples
        self.dataMap = dataMap
        self.gPathway = gPathway
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        
        ## copy tables
        system("echo Output Reports ... >> progress.log")
        system("cp param_%s/sig.tab sig.tab" % ("_".join(self.topAUC_params)))
        system("cp param_%s/pshift.train.tab pshift.tab" % ("_".join(self.topAUC_params)))
        
        ## output msep and background plots
        system("%s/msep.R %s param_%s/mut.train.scores param_%s/non.train.scores" % (self.directory, 
                                                 self.mutatedGene, "_".join(self.topAUC_params), 
                                                 "_".join(self.topAUC_params)))
        system("%s/background.R %s param_%s/real.scores param_%s/null.scores" % (self.directory, 
                                                 self.mutatedGene, "_".join(self.topAUC_params), 
                                                 "_".join(self.topAUC_params)))
        
        ## prepare web report
        system("mkdir img")
        
        ## draw circles
        system("cat up_%s | transpose.pl > param_%s/exp.tab" % (re.split("/", self.dataMap["exp"])[-1], "_".join(self.topAUC_params)))
        os.chdir("param_%s" % ("_".join(self.topAUC_params)))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f mut.features ../img/ mut.circle exp.tab paradigm_up.tab paradigm_down.tab shift.circle" % (circleExec, self.mutatedGene))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f up.features ../img/ mut.circle exp.tab paradigm_up.tab" % (circleExec, self.mutatedGene))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f down.features ../img/ mut.circle exp.tab paradigm_down.tab" % (circleExec, self.mutatedGene))
        os.chdir("..")
        
        ## output sif
        (uNodes, uInteractions) = rPathway("param_%s/upstream_pathway.tab" % 
                                                    ("_".join(self.topAUC_params)))
        (dNodes, dInteractions) = rPathway("param_%s/downstream_pathway.tab" % 
                                                    ("_".join(self.topAUC_params)))
        cPathway = combinePathways(Pathway(uNodes, uInteractions), Pathway(dNodes, dInteractions))
        wSIF("pshift_%s.sif" % (self.mutatedGene), cPathway.interactions)
        
        ## node attributes
        scoreMap = {}
        scoreMap["pshift_%s" % (self.mutatedGene)] = {}
        for node in self.gPathway.nodes.keys():
            if node == self.mutatedGene:
                scoreMap["pshift_%s" % (self.mutatedGene)][node] = 10
            else:
                scoreMap["pshift_%s" % (self.mutatedGene)][node] = 7
        wNodeAttributes(self.gPathway.nodes, scoreMap = scoreMap, directory = "./")

class pshiftReport(Target):
    def __init__(self, includeFeatures, reportDir, directory):
        Target.__init__(self, time=10000)
        self.includeFeatures = includeFeatures
        self.reportDir = reportDir
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## cytoscape-web
        for mutatedGene in self.includeFeatures:
            if os.path.exists("analysis/%s/sig.tab" % (mutatedGene)):
                tableFiles = []
                tableFiles.append("analysis/%s/sig.tab" % (mutatedGene))
                tableFiles.append("msepPlot:analysis/%s/%s.msep.pdf" % (mutatedGene, mutatedGene))
                tableFiles.append("backgroundPlot:analysis/%s/%s.background.pdf" % (mutatedGene, mutatedGene))
                tableFiles.append("analysis/%s/avgAUC.tab" % (mutatedGene))
                tableFiles.append("analysis/%s/pshift.tab" % (mutatedGene))
                system("pathmark-report.py -t %s analysis/%s %s" % (",".join(tableFiles), mutatedGene, self.reportDir))
                system("cp analysis/%s/pshift* %s" % (mutatedGene, self.reportDir))

def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory detected, remove it first to start a fresh run"
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] paradigmDir mutFile ")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help = "Add as child of jobFile rather than a new jobTree")
    parser.add_option("-s", "--samples", dest="includeSamples", default="")
    parser.add_option("-f", "--features", dest="includeFeatures", default="")
    options, args = parser.parse_args()
    print "Using Batch System '%s'" % (options.batchSystem)
    
    if len(args) == 1:
        if args[0] == "clean":
            cmd = "rm -rf .jobTree analysis *.R"
            print cmd
            os.system(cmd)
            sys.exit(0)
    
    assert len(args) == 2
    paradigmDir = os.path.abspath(args[0])
    mutFile = args[1]
    sampleFile = options.includeSamples
    featureFile = options.includeFeatures
    
    ## paramMap
    paramMap = {}
    paramMap["dist"] = [2]
    paramMap["thresh"] = [0.68]
    paramMap["inc"] = [0.01]
    paramMap["method"] = ["vsNon"]
    paramMap["stat"] = "tt"
    if os.path.exists("mut.cfg"):
        f = open("mut.cfg", "r")
        for line in f:
            if line.isspace():
                continue
            pline = re.split("\t", line.rstrip("\r\n"))
            if line.startswith("distanceParams"):
                paramMap["dist"] = [int(i) for i in re.split(",", pline[1])]
            elif line.startswith("tstatParams"):
                paramMap["thresh"] = [float(i) for i in re.split(",", pline[1])]
            elif line.startswith("incrParams"):
                paramMap["inc"] = [float(i) for i in re.split(",", pline[1])]
            elif line.startswith("methodParams"):
                paramMap["method"] = [i for i in re.split(",", pline[1])]
            elif line.startswith("signalMethod"):
               paramMap["stat"] = int(pline[1])
            elif line.startswith("cohortName"):
                paramMap["cohortName"] = pline[1]
        f.close()
    if "cohortName" not in paramMap:
        paramMap["cohortName"] = re.split("/", os.getcwd())[-1]
    
    ## foldMap
    rRepeats = 1
    mFolds = 5
    foldMap = {}
    if os.path.exists("fold.map"):
        (mapData, mapSamples, mapRepeats) = rCRSData("fold.map", retFeatures = True)
        for r, row in enumerate(mapRepeats):
            foldMap[r+1] = {}
            for col in mapSamples:
                mFolds = max(mFolds, int(mapData[col][row]))
            for m in range(mFolds):
                foldMap[r+1][m+1] = []
            for col in mapSamples:
                foldMap[r+1][int(mapData[col][row])].append(col)
            if (r+1) == rRepeats:
                break
    else:
        for r in range(1, rRepeats+1):
            foldMap[r] = {}
            for m in range(1, mFolds+1):
                foldMap[r][m] = None
    
    ## check files
    pathwayFile = None
    cnvFile = None
    expFile = None
    dataMap = {}
    assert os.path.exists("%s/clusterFiles" % (paradigmDir))
    for file in os.listdir("%s/clusterFiles" % (paradigmDir)):
        if file.endswith("pathway.tab"):
            pathwayFile = "%s/clusterFiles/%s" % (paradigmDir, file)
        elif file.endswith("CNV.tab"):
            cnvFile = "%s/clusterFiles/%s" % (paradigmDir, file)
            dataMap["cnv"] = cnvFile
        elif (file.endswith("Expression.tab") | file.endswith("Expression.vCohort.tab") | 
              file.endswith("Expression.vNormal.tab")):
            expFile = "%s/clusterFiles/%s" % (paradigmDir, file)
            dataMap["exp"] = expFile
    assert (pathwayFile != None)
    assert (cnvFile != None)
    assert (expFile != None)
    paradigmFile = None
    if os.path.exists("%s/merge_merged_unfiltered.tab" % (paradigmDir)):
        paradigmFile = "%s/merge_merged_unfiltered.tab" % (paradigmDir)
    elif os.path.exists("%s/merge_merged.tab" % (paradigmDir)):
        paradigmFile = "%s/merge_merged.tab" % (paradigmDir)
	
    ## store feature, sample and pathway information
    dataFeatures = list(set(retColumns(cnvFile)) & set(retColumns(expFile)))
    includeFeatures = None
    if len(featureFile) != 0:
        includeFeatures = rList(featureFile)
    
    dataSamples = list(set(retRows(cnvFile)) & set(retRows(expFile)))
    if len(sampleFile) != 0:
        dataSamples = list(set(dataSamples) & set(rList(sampleFile)))

    (gNodes, gInteractions) = rPathway(pathwayFile)
    gfPathway = flattenPathway(Pathway(gNodes, gInteractions))
    if not useFlattened:
        gPathway = Pathway(gNodes, gInteractions)
    else:
        gPathway = gfPathway
	
    mutationOrder = []
    mutationMap = {}
    f = open(mutFile, "r")
    for line in f:
        if line.isspace():
            continue
        pline = re.split("\t", line.rstrip("\r\n"))
        mutatedGene = pline[0]
        mutatedSamples = list(set(re.split(",", pline[2])) & set(dataSamples))
        if mutatedGene in gPathway.nodes:
            if len(mutatedSamples) >= mutationThreshold:
                mutationMap[mutatedGene] = deepcopy(mutatedSamples)
                if includeFeatures is None:
                    mutationOrder.append(mutatedGene)
    f.close()
    if includeFeatures is not None:
        for mutatedGene in includeFeatures:
            if mutatedGene in mutationMap:
                mutationOrder.append(mutatedGene)
     
    submitMap = {}
    for mutatedGene in mutationOrder:
        submitMap[mutatedGene] = deepcopy(mutationMap[mutatedGene])
        if len(submitMap.keys()) >= maxFeatures:
            break
    
    ## run
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    
    s = Stack(branchGenes(dataSamples, dataFeatures, dataMap, submitMap, gPathway, paradigmDir, 
              paramMap, foldMap, os.getcwd()))
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
    main()     
