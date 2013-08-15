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

## select seed
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
runParallel = False              ## Run different genes in parallel
mutationThreshold = 5            ## Minimum mutations in cohort needed to consider a gene
maxFeatures = 30                 ## Maximum number of features to in PARADIGM Shift run

## object classes
class ParadigmSetup:
    def __init__(self, directory, include, public = False, size = 50, nulls = 30):
        self.directory = directory
        self.config = None
        self.params = None
        self.pathway = None
        self.cnv = None
        self.exp = None
        self.prot = None
        self.act = None
        self.ipl = None
        self.features = []
        self.samples = []
        assert os.path.exists("%s/config.txt" % (directory))
        self.config = "%s/config.txt" % (directory)
        assert os.path.exists("%s/params.txt" % (directory))
        self.params = "%s/params.txt" % (directory)        
        assert os.path.exists("%s/clusterFiles" % (directory))
        for file in os.listdir("%s/clusterFiles" % (directory)):
            if file.endswith("pathway.tab"):
                self.pathway = "%s/clusterFiles/%s" % (directory, file)
        assert(self.pathway != None)
        f = open(self.config, "r")
        for line in f:
            if line.isspace():
                continue
            if line.startswith("evidence"):
                tokens = line.lstrip("evidence [").rstrip("]").split(",")
                attachNode = None
                attachFile = None
                for token in tokens:
                    parts = token.split("=")
                    if parts[0] == "node":
                        attachNode = parts[1]
                    elif parts[0] == "suffix":
                        attachFile = parts[1]
                if attachNode == "genome":
                    self.cnv = "%s/clusterFiles/%s" % (directory, attachFile)
                elif attachNode == "mRNA":
                    self.exp = "%s/clusterFiles/%s" % (directory, attachFile)
                elif attachNode == "protein":
                    self.prot = "%s/clusterFiles/%s" % (directory, attachFile)
                elif attachNode == "active":
                    self.act = "%s/clusterFiles/%s" % (directory, attachFile)
        f.close()
        assert(self.exp != None)
        if os.path.exists("%s/merge_merged_unfiltered.tab" % (directory)):
            self.ipl = "%s/merge_merged_unfiltered.tab" % (directory)
        elif os.path.exists("%s/merge_merged.tab" % (directory)):
            self.ipl = "%s/merge_merged.tab" % (directory)
        self.features = list(set(retColumns(self.cnv)) & set(retColumns(self.exp)))
        self.samples = list(set(retRows(self.cnv)) & set(retRows(self.exp)))
        if include:
            self.samples = list(set(self.samples) & set(rList(include)))
        self.public = public
        self.size = size
        self.nulls = nulls
        
class Parameters:
    def __init__(self):
        self.distance = [2]
        self.threshold = [0.68]
        self.penalty = [0.01]
        self.selection = ["vsNon"]
        self.override = [None, None]
        self.cv = True
        self.msep = "tt"
        self.cohort = os.getcwd().split("/")[-1]
        self.link = False
        self.clean = True
    def set(self, configFile):
        f = open(configFile, "r")
        for line in f:
            if line.isspace():
                continue
            pline = line.rstrip().split("\t")
            if line.startswith("distance"):
                self.distance = [int(item) for item in pline[1].split(",")]
            elif line.startswith("threshold"):
                self.threshold = [float(item) for item in pline[1].split(",")]
            elif line.startswith("penalty"):
                self.penalty = [float(item) for item in pline[1].split(",")]
            elif line.startswith("selection"):
                self.selection = [item for item in pline[1].split(",")]
            elif line.startswith("override"):
                self.override = []
                for item in pline[1].split(","):
                    if len(item) == 0:
                        self.override.append(None)
                    else:
                        self.override.append(os.path.abspath(item))
            elif line.startswith("cv"):
                self.cv = bool(int(pline[1]))
            elif line.startswith("msep"):
                self.msep = pline[1]
            elif line.startswith("cohort"):
                self.cohort = pline[1]
            elif line.startswith("search"):
                self.link = bool(int(pline[1]))
            elif line.startswith("clean"):
                self.clean = bool(int(pline[1]))
        f.close()

class Mutated:
    def __init__(self, gene, all, positive, negative = None, directory = ""):
        self.gene = gene
        self.positive = list(set(positive) & set(all))
        if negative:
            self.negative = list(set(negative) & set(all))
        else:
            self.negative = list(set(all) - set(positive))
        self.data = "%s/analysis/%s/data/" % (directory, gene)
    def proportion(self):
        return(float(len(self.positive))/float(len(self.positive) + len(self.negative)))

def batchCount(cohortSize, batchSize = 15):
    return(int((cohortSize+batchSize-1)/batchSize))

def chunks(sampleList, batchSize = 15):
    for i in xrange(0, len(sampleList), batchSize):
        yield sampleList[i:i+batchSize]

## R script definitions
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

## jobTree classes
class jtData(Target):
    def __init__(self, mutatedGene, useSamples, useFeatures, allFeatures, dataFiles,
                 directory, paradigmPublic = False, batchSize = 15):
        Target.__init__(self, time=1000)
        self.mutatedGene = mutatedGene
        self.useSamples = useSamples
        self.useFeatures = useFeatures
        self.allFeatures = allFeatures
        self.dataFiles = dataFiles
        self.directory = directory
        self.paradigmPublic = paradigmPublic
        self.batchSize = batchSize
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        random.seed(randSeed)
        
        ## write data
        for file in self.dataFiles:
            rwCRSData("data/full_%s" % (file.split("/")[-1]), file, useRows = self.useSamples)
            rwCRSData("data/up_%s" % (file.split("/")[-1]), file, useRows = self.useSamples, useCols = self.useFeatures)
            rwCRSData("data/down_%s" % (file.split("/")[-1]), file, useRows = self.useSamples, useCols = list(set(self.useFeatures)-set([self.mutatedGene])))
        
        ## public batching
        if self.paradigmPublic:
            batches = batchCount(len(self.useSamples), batchSize = self.batchSize)
            batchChunker = chunks(self.useSamples, batchSize = self.batchSize)
            for b in range(batches):
                currentSamples = batchChunker.next()
                for file in self.dataFiles:
                    rwCRSData("data/up_b%s_%s_%s" % (b, batches, file.split("/")[-1]), file, useRows = currentSamples, useCols = self.useFeatures)
                    rwCRSData("data/down_b%s_%s_%s" % (b, batches, file.split("/")[-1]), file, useRows = currentSamples, useCols = list(set(self.useFeatures)-set([self.mutatedGene])))

class jtNData(Target):
    def __init__(self, null, mutatedGene, useSamples, useFeatures, allFeatures, dataFiles,
                 directory, paradigmPublic = False, batchSize = 15):
        Target.__init__(self, time=1000)
        self.null = null
        self.mutatedGene = mutatedGene
        self.useSamples = useSamples
        self.useFeatures = useFeatures
        self.allFeatures = allFeatures
        self.dataFiles = dataFiles
        self.directory = directory
        self.paradigmPublic = paradigmPublic
        self.batchSize = batchSize
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedGene)
        os.chdir(shiftDir)
        random.seed(randSeed+self.null)
        
        ## permute columns for nulls
        colMap = {}
        colFeatures = list(set(self.allFeatures)-set([self.mutatedGene]))
        colPermute = random.sample(colFeatures, len(colFeatures))
        for i in range(0, len(colFeatures)):
            colMap[colFeatures[i]] = colPermute[i]
        
        ## write data
        for file in self.dataFiles:
            rwCRSData("data/up_N%s_%s" % (self.null, file.split("/")[-1]), file, useRows = self.useSamples, colMap = colMap, useCols = self.useFeatures)
            rwCRSData("data/down_N%s_%s" % (self.null, file.split("/")[-1]), file, useRows = self.useSamples, colMap = colMap, useCols = list(set(self.useFeatures)-set([self.mutatedGene])))

        ## public batching
        if self.paradigmPublic:
            batches = batchCount(len(self.useSamples), batchSize = self.batchSize)
            batchChunker = chunks(self.useSamples, batchSize = self.batchSize)
            for b in range(batches):
                currentSamples = batchChunker.next()
                for file in self.dataFiles:
                    rwCRSData("data/up_N%s_b%s_%s_%s" % (self.null, b, batches, file.split("/")[-1]), file, useRows = currentSamples, colMap = colMap, useCols = self.useFeatures)
                    rwCRSData("data/down_N%s_b%s_%s_%s" % (self.null, b, batches, file.split("/")[-1]), file, useRows = currentSamples, colMap = colMap, useCols = list(set(self.useFeatures)-set([self.mutatedGene])))

class jtCmd(Target):
    def __init__(self, cmd, directory):
        Target.__init__(self, time=1000)
        self.cmd = cmd
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        system(self.cmd)

class queueGenes(Target):
    def __init__(self, mutatedList, paradigmSetup, gPathway, params, foldMap, directory,
                 appendFeatures = []):
        Target.__init__(self, time=10000)
        self.mutatedList = mutatedList
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.params = params
        self.foldMap = foldMap
        self.directory = directory
        self.appendFeatures = appendFeatures
    def run(self):
        os.chdir(self.directory)
        
        ## queue genes
        htmlFeatures = deepcopy(self.appendFeatures)
        if not os.path.exists("analysis"):
            system("mkdir analysis")
        if len(self.mutatedList) > 0:
            mutatedFeature = self.mutatedList[0]
            if not os.path.exists("analysis/%s" % (mutatedFeature.gene)):
                system("mkdir analysis/%s" % (mutatedFeature.gene))
                htmlFeatures.append(mutatedFeature.gene)
                self.addChildTarget(branchFolds(mutatedFeature, self.paradigmSetup,
                                                self.gPathway, self.params,
                                                self.foldMap, self.directory))
            self.setFollowOnTarget(queueGenes(self.mutatedList[1:], self.paradigmSetup,
                                              self.gPathway, self.params, self.foldMap,
                                              self.directory,
                                              appendFeatures = htmlFeatures))
        else:
            if os.path.exists(htmlDir):
                self.setFollowOnTarget(pshiftReport(htmlFeatures, "%s/%s" % (htmlDir, self.params.cohort), self.directory))

class branchGenes(Target):
    def __init__(self, mutatedList, paradigmSetup, gPathway, params, foldMap, directory):
        Target.__init__(self, time=10000)
        self.mutatedList = mutatedList
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.params = params
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## branch genes
        htmlFeatures = []
        if not os.path.exists("analysis"):
            system("mkdir analysis")
        for mutatedFeature in self.mutatedList:
            if not os.path.exists("analysis/%s" % (mutatedFeature.gene)):
                system("mkdir analysis/%s" % (mutatedFeature.gene))
                htmlFeatures.append(mutatedFeature.gene)
                self.addChildTarget(branchFolds(mutatedFeature, self.paradigmSetup,
                                                self.gPathway, self.params,
                                                self.foldMap, self.directory))
        if os.path.exists(htmlDir):
            self.setFollowOnTarget(pshiftReport(htmlFeatures, "%s/%s" % (htmlDir, self.params.cohort), self.directory))

class branchFolds(Target):
    def __init__(self, mutatedFeature, paradigmSetup, gPathway, params, foldMap, directory):
        Target.__init__(self, time=10000)
        self.mutatedFeature = mutatedFeature
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.params = params
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedFeature.gene)
        os.chdir(shiftDir)
        random.seed(randSeed+123454321)
        
        ## prepare htmlDir
        
        ## get proteinFeatures based on maxDist and gPathway
        system("echo Preparing Genomic Data ... >> progress.log")
        proteinFeatures = []
        for feature in getNeighbors(self.mutatedFeature.gene, max(self.params.distance) + 1, self.gPathway.interactions):
            if self.gPathway.nodes[feature] == "protein":
                proteinFeatures.append(feature)
        
        ## write data
        system("mkdir data")
        dataFiles = []
        if self.paradigmSetup.cnv is not None:
            dataFiles.append(self.paradigmSetup.cnv)
        if self.paradigmSetup.exp is not None:
            dataFiles.append(self.paradigmSetup.exp)
        if self.paradigmSetup.prot is not None:
            dataFiles.append(self.paradigmSetup.prot)
        if self.paradigmSetup.act is not None:
            dataFiles.append(self.paradigmSetup.act)
        genData = jtData(self.mutatedFeature.gene, self.paradigmSetup.samples, proteinFeatures, self.paradigmSetup.features, 
                         dataFiles, self.directory, paradigmPublic = self.paradigmSetup.public,
                         batchSize = self.paradigmSetup.size)
        genData.run()
        for null in range(1, self.paradigmSetup.nulls+1):
            genData = jtNData(null, self.mutatedFeature.gene, self.paradigmSetup.samples, proteinFeatures, self.paradigmSetup.features,
                              dataFiles, self.directory, paradigmPublic = self.paradigmSetup.public,
                              batchSize = self.paradigmSetup.size)
            genData.run()
        
        ## pick samples
        system("echo Branching Folds and Params... >> progress.log")
        rRepeats = len(self.foldMap.keys())
        mFolds = len(self.foldMap[1].keys())
        if self.foldMap[1][1] is None:
            mutSamples = self.mutatedFeature.positive
            nonSamples = self.mutatedFeature.negative
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
        if self.params.cv:
            for r in range(1, rRepeats+1):
                for f in range(1, mFolds+1):
                    fold = (r-1)*mFolds+f
                    system("mkdir fold%s" % (fold))
                    self.addChildTarget(branchParams(fold, self.mutatedFeature,
                                                     foldSamples[r][f], self.paradigmSetup,
                                                     self.gPathway, self.params,
                                                     self.foldMap, self.directory))
        
        ## run final
        self.setFollowOnTarget(branchParams(0, self.mutatedFeature, 
                                            self.paradigmSetup.samples, self.paradigmSetup,
                                            self.gPathway, self.params, self.foldMap,
                                            self.directory))

class branchParams(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, gPathway,
                 params, foldMap, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.params = params
        self.foldMap = foldMap
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedFeature.gene)
            os.chdir(shiftDir)
            system("echo Building Final Model ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutatedFeature.gene, self.fold)
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
                    if os.path.exists("fold%s/auc.stat" % (fold)):
                        f = open("fold%s/auc.stat" % (fold), "r")
                        line = f.readline()
                        f.close()
                        (auc_tr, auc_te, mparams) = line.rstrip().split("\t")
                    else:
                        (auc_tr, auc_te, mparams) = ("---", "---", "---")
                    aucList.append(auc_te)
                    aucLines.append("%s\t%s\t%s\t%s\n" % (fold, auc_tr, auc_te, mparams))
            o = open("avgAUC.tab", "w")
            o.write("> %s\tavgAUC:%s\n" % (self.mutatedFeature.gene, mean(aucList)))
            o.write("# fold\ttrain\ttest\tparams\n")
            for line in aucLines:
                o.write(line)
            o.close()
        
        ## branch params
        for distance in self.params.distance:
            for thresh in self.params.threshold:
                for inc in self.params.penalty:
                    for method in self.params.selection:
                        system("mkdir param_%s_%s_%s_%s" % (distance, thresh, inc, method))
                        os.chdir("param_%s_%s_%s_%s" % (distance, thresh, inc, method))
                        system("cp %s/config.txt config.txt" % (self.paradigmSetup.directory))
                        system("cp %s/params.txt params.txt" % (self.paradigmSetup.directory))
                        for file in os.listdir(self.paradigmSetup.directory):
                            if file.endswith(".imap"):
                                system("cp %s/%s %s" % (self.paradigmSetup.directory, file, file))
                            elif file.endswith(".dogma"):
                                system("cp %s/%s %s" % (self.paradigmSetup.directory, file, file))
                        os.chdir("..")
                        self.addChildTarget(prepareNeighborhood(self.fold, 
                                        self.mutatedFeature, self.trainSamples,
                                        self.paradigmSetup, distance, thresh, inc, method,
                                        self.gPathway, self.directory,
                                        setUpstream = self.params.override[0],
                                        setDownstream = self.params.override[1],
                                        linkSearch = self.params.link))
        
        ## evaluate models
        self.setFollowOnTarget(evaluateParams(self.fold, self.mutatedFeature, self.trainSamples, 
                                              self.paradigmSetup, self.gPathway, 
                                              self.params, self.directory))
        
class prepareNeighborhood(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, dParam,
                 tParam, iParam, method, gPathway, directory,
                 setUpstream = None, setDownstream = None, linkSearch = False):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.dParam = dParam
        self.tParam = tParam
        self.iParam = iParam
        self.method = method
        self.gPathway = gPathway
        self.directory = directory
        self.setUpstream = setUpstream
        self.setDownstream = setDownstream
        self.linkSearch = linkSearch
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s/param_%s_%s_%s_%s" % (self.directory, self.mutatedFeature.gene, 
                                                             self.dParam, self.tParam, self.iParam, 
                                                             self.method)
            os.chdir(shiftDir)
            system("echo Preparing Network ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s/param_%s_%s_%s_%s" % (self.directory,
                                                                    self.mutatedFeature.gene, self.fold, 
                                                                    self.dParam, self.tParam, 
                                                                    self.iParam, self.method)
            os.chdir(shiftDir)
        
        ## get upstream and downstream
        system("echo Selecting Mutational Network ... >> progress.log")
        dataFile = "%s/full_%s" % (self.mutatedFeature.data, self.paradigmSetup.exp.split("/")[-1])
        positiveSamples = list(set(self.trainSamples) & set(self.mutatedFeature.positive))
        negativeSamples = list(set(self.trainSamples) & set(self.mutatedFeature.negative))
        (featureScore, featureRank) = getFeatureScores(dataFile, positiveSamples, 
                                               negativeSamples, method = self.method,
                                               outFile = "feature.score")
        focusPathways = identifyNeighborhoods(self.mutatedFeature.gene, self.gPathway,
                                              maxDistance = self.dParam)
        focusScore = {}
        for index, focusPathway in enumerate(focusPathways):
            upstreamFeatures = focusPathway[0][0]
            downstreamFeatures = focusPathway[1][0]
            featureValues = []
            for feature in upstreamFeatures:
                if feature in featureScore:
                    featurePaths = shortestPath(feature, self.mutatedFeature.gene, focusPathway[0][1].interactions)
                    pathSigns = [signPath(path, focusPathway[0][1].interactions) for path in featurePaths]
                    try:
                        featureValues.append((-1)*featureScore[feature]*mean(pathSigns))
                    except TypeError:
                        pass
            for feature in downstreamFeatures:
                if feature in featureScore:
                    featurePaths = shortestPath(self.mutatedFeature.gene, feature, focusPathway[1][1].interactions)
                    pathSigns = [signPath(path, focusPathway[1][1].interactions) for path in featurePaths]
                    try:
                        featureValues.append((1)*featureScore[feature]*mean(pathSigns))
                    except TypeError:
                        pass
            if len(featureValues) > 0:
                focusScore[index] = mean(featureValues)
            else:
                focusScore[index] = 0.0
        rankFocus = list(set(focusScore.keys()) - set([0]))
        rankFocus.sort(lambda x, y: cmp(focusScore[x],focusScore[y]))
        addedFocus = [0]
        #### while (float(len(addedFocus))/float(len(focusScore.keys())) < 0.5) or (len(addedFocus) < 4):
        while len(rankFocus) > 0:
            addedFocus.append(rankFocus.pop(0))
        uFeatures = []
        dFeatures = []
        uBase = Pathway({self.mutatedFeature.gene : self.gPathway.nodes[self.mutatedFeature.gene]}, {})
        dBase = Pathway({self.mutatedFeature.gene : self.gPathway.nodes[self.mutatedFeature.gene]}, {})
        for focus in addedFocus:
            system("echo ... %s >> progress.log" % (focusScore[focus]))
            uFeatures = list(set(uFeatures) | set(focusPathways[focus][0][0]))
            dFeatures = list(set(dFeatures) | set(focusPathways[focus][1][0]))
            uBase = combinePathways(uBase, focusPathways[focus][0][1])
            dBase = combinePathways(dBase, focusPathways[focus][1][1])
        uFeatures = list(set(uFeatures) - set(dFeatures))
        wPathway("upstream_base.tab", uBase.nodes, uBase.interactions)
        wPathway("downstream_base.tab", dBase.nodes, dBase.interactions)
        if len(uFeatures) >= 4 and (dFeatures) >= 4:
            isPass = True
            (uPathway, dPathway) = selectMutNeighborhood(self.mutatedFeature.gene, featureRank, uFeatures, dFeatures, uBase, dBase, threshold = self.tParam, penalty = self.iParam)
            if self.setUpstream is not None:
                (uNodes, uInteractions) = rPathway(self.setUpstream)
                uPathway = Pathway(uNodes, uInteractions)
            if self.setDownstream is not None:
                (dNodes, dInteractions) = rPathway(self.setDownstream)
                dPathway = Pathway(dNodes, dInteractions)
            wPathway("upstream_pathway.tab", uPathway.nodes, uPathway.interactions)
            wPathway("downstream_pathway.tab", dPathway.nodes, dPathway.interactions)
            
            uSelected = list(set(uFeatures) & set(uPathway.nodes))
            dSelected = list(set(dFeatures) & set(dPathway.nodes))
        else:
            isPass = False
        
        ## check thresholds
        if not isPass:
            system("echo Stop ... >> progress.log")
            f = open("auc.stat", "w")
            f.write("---\t---\n")
            f.close()
        
        elif self.linkSearch:
            self.setFollowOnTarget(runLinkSearch(self.fold, self.mutatedFeature,
                                                 self.trainSamples, self.paradigmSetup,
                                                 featureRank, self.tParam, self.iParam,
                                                 uPathway, dPathway,
                                                 uBase, dBase, self.directory, shiftDir))
        else:
            self.setFollowOnTarget(runPARADIGM(self.fold, self.mutatedFeature,
                                                self.trainSamples, self.paradigmSetup,
                                                uPathway, dPathway, shiftDir))

class runLinkSearch(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, featureRank,
                 tParam, iParam, uPathway, dPathway, uBase, dBase, root, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.featureRank = featureRank
        self.tParam = tParam
        self.iParam = iParam
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.uBase = uBase
        self.dBase = dBase
        self.root = root
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## assert links.tab exists
        assert os.path.exists("%s/links.tab" % (self.root))
        (lNodes, lInteractions) = rPathway("%s/links.tab" % (self.root))
        
        ## determine possible additions
        upstreamAdditions = []
        downstreamAdditions = []
        for source in lInteractions:
            for target in lInteractions[source]:
                if target in self.uBase.nodes:
                    if lInteractions[source][target].startswith("-a"):
                        if source in self.uBase.interactions:
                            if target not in self.uBase.interactions[source]:
                                upstreamAdditions.append( (source, target, lInteractions[source][target]) )
                        elif source not in self.uBase.nodes:
                            upstreamAdditions.append( (source, target, lInteractions[source][target]) )
                if source in self.dBase.nodes:
                    if lInteractions[source][target].startswith("-t"):
                        if source in self.dBase.interactions:
                            if target not in self.dBase.interactions[source]:
                                downstreamAdditions.append( (source, target, lInteractions[source][target]) )
                        else:
                            downstreamAdditions.append( (source, target, lInteractions[source][target]) )
        ## score and determine tested additions
        additionMap = {}
        additionRank = {}
        index = 1
        for addition in upstreamAdditions:
            log("Candidate: %s, %s, %s\n" % addition, file = "addition.log")
            if addition[0] not in self.featureRank:
                continue
            elif self.featureRank[addition[0]] >= self.tParam:
                log("u%s - %s\n" % (index, addition), file = "addition.log")
                additionMap["u%s" % (index)] = addition
                additionRank["u%s" % (index)] = self.featureRank[addition[0]]
                index += 1
        index = 1
        for addition in downstreamAdditions:
            log("Candidate: %s, %s, %s\n" % addition, file = "addition.log")
            if addition[1] not in self.featureRank:
                continue
            elif self.featureRank[addition[1]] >= self.tParam:
                log("d%s - %s\n" % (index, addition), file = "addition.log")
                additionMap["d%s" % (index)] = addition
                additionRank["d%s" % (index)] = self.featureRank[addition[1]]
                index += 1
        rankAddition = additionMap.keys()
        rankAddition.sort(lambda x, y: cmp(additionRank[y],additionRank[x]))
        if len(rankAddition) > 15:
            for i in rankAddition[15:]:
                if i in additionMap:
                    del additionMap[i]
        
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
        self.addChildTarget(runPARADIGM(self.fold, self.mutatedFeature, 
                                        self.trainSamples, self.paradigmSetup, self.uPathway, 
                                        self.dPathway, "%s/iter_0" % (self.directory)))
        self.setFollowOnTarget(branchLinkSearch(self.fold, self.mutatedFeature,
                                        self.trainSamples, self.paradigmSetup,
                                        self.uPathway, self.dPathway,
                                        self.uBase, self.dBase, Pathway(lNodes, lInteractions),
                                        1, 5, [], additionMap, self.directory))

class branchLinkSearch(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, uPathway, 
                 dPathway, uBase, dBase, extensions, iter, maxiter, additions, additionMap, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.uBase = uBase
        self.dBase = dBase
        self.extensions = extensions
        self.iter = iter
        self.maxiter = maxiter
        self.additions = additions
        self.additionMap = additionMap
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## get current pathway
        currentUp = deepcopy(self.uPathway)
        currentDown = deepcopy(self.dPathway)
        for addition in self.additions:
            if addition.startswith("u"):
                if self.additionMap[addition][1] not in currentUp.nodes:
                    addPaths = shortestPath(self.additionMap[addition][1], self.mutatedFeature.gene, self.uBase.interactions)
                    for path in addPaths:
                        for i in range(len(path)-1):
                            currentUp.nodes[path[i]] = self.uBase.nodes[path[i]]
                            if path[i] not in currentUp.interactions:
                                currentUp.interactions[path[i]] = {}
                            currentUp.interactions[path[i]][path[i+1]] = self.uBase.interactions[path[i]][path[i+1]]
                if self.additionMap[addition][0] not in currentUp.nodes:
                    currentUp.nodes[self.additionMap[addition][0]] = self.extensions.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][0] not in currentUp.interactions:
                    currentUp.interactions[self.additionMap[addition][0]] = {}
                currentUp.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            else:
                if self.additionMap[addition][0] not in currentDown.nodes:
                    addPaths = shortestPath(self.mutatedFeature.gene, self.additionMap[addition][0], self.dBase.interactions)
                    for path in addPaths:
                        for i in range(len(path)-1):
                            currentDown.nodes[path[i+1]] = self.dBase.nodes[path[i+1]]
                            if path[i] not in currentDown.interactions:
                                currentDown.interactions[path[i]] = {}
                            currentDown.interactions[path[i]][path[i+1]] = self.dBase.interactions[path[i]][path[i+1]]
                if self.additionMap[addition][1] not in currentDown.nodes:
                    currentDown.nodes[self.additionMap[addition][1]] = self.extensions.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in currentDown.interactions:
                    currentDown.interactions[self.additionMap[addition][0]] = {}
                currentDown.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
        
        ## for each addition
        for addition in self.additionMap:
            if addition in self.additions:
                continue
            addUp = deepcopy(currentUp)
            addDown = deepcopy(currentDown)
            if addition.startswith("u"):
                if self.additionMap[addition][1] not in addUp.nodes:
                    addPaths = shortestPath(self.additionMap[addition][1], self.mutatedFeature.gene, self.uBase.interactions)
                    for path in addPaths:
                        for i in range(len(path)-1):
                            addUp.nodes[path[i]] = self.uBase.nodes[path[i]]
                            if path[i] not in addUp.interactions:
                                addUp.interactions[path[i]] = {}
                            addUp.interactions[path[i]][path[i+1]] = self.uBase.interactions[path[i]][path[i+1]]
                if self.additionMap[addition][0] not in addUp.nodes:
                    addUp.nodes[self.additionMap[addition][0]] = self.extensions.nodes[self.additionMap[addition][0]]
                if self.additionMap[addition][0] not in addUp.interactions:
                    addUp.interactions[self.additionMap[addition][0]] = {}
                addUp.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            else:
                if self.additionMap[addition][0] not in addDown.nodes:
                    addPaths = shortestPath(self.mutatedFeature.gene, self.additionMap[addition][0], self.dBase.interactions)
                    for path in addPaths:
                        for i in range(len(path)-1):
                            addDown.nodes[path[i+1]] = self.dBase.nodes[path[i+1]]
                            if path[i] not in addDown.interactions:
                                addDown.interactions[path[i]] = {}
                            addDown.interactions[path[i]][path[i+1]] = self.dBase.interactions[path[i]][path[i+1]]
                if self.additionMap[addition][1] not in addDown.nodes:
                    addDown.nodes[self.additionMap[addition][1]] = self.extensions.nodes[self.additionMap[addition][1]]
                if self.additionMap[addition][0] not in addDown.interactions:
                    addDown.interactions[self.additionMap[addition][0]] = {}
                addDown.interactions[self.additionMap[addition][0]][self.additionMap[addition][1]] = self.additionMap[addition][2]
            system("mkdir iter_%s-%s" % (self.iter, addition))
            wPathway("iter_%s-%s/upstream_pathway.tab" % (self.iter, addition), addUp.nodes, addUp.interactions)
            wPathway("iter_%s-%s/downstream_pathway.tab" % (self.iter, addition), addDown.nodes, addDown.interactions)
            system("cp config.txt iter_%s-%s/config.txt" % (self.iter, addition))
            system("cp params.txt iter_%s-%s/params.txt" % (self.iter, addition))
            for file in os.listdir("."):
                if file.endswith(".imap"):
                    system("cp %s iter_%s-%s/%s" % (file, self.iter, addition, file))
                elif file.endswith(".dogma"):
                    system("cp %s iter_%s-%s/%s" % (file, self.iter, addition, file))
            self.addChildTarget(runPARADIGM(self.fold, self.mutatedFeature, 
                                            self.trainSamples, self.paradigmSetup, addUp, 
                                            addDown, 
                                            "%s/iter_%s-%s" % (self.directory, self.iter, addition)))
        self.setFollowOnTarget(selectLinkOptimum(self.fold, self.mutatedFeature, 
                                                 self.trainSamples, self.paradigmSetup,
                                                 self.uPathway, self.dPathway, self.uBase,
                                                 self.dBase, self.extensions, self.iter,
                                                 self.maxiter, self.additions,
                                                 self.additionMap, self.directory))

class selectLinkOptimum(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, uPathway,
                 dPathway, uBase, dBase, extensions, iter, maxiter, additions, additionMap, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.uBase = uBase
        self.dBase = dBase
        self.extensions = extensions
        self.iter = iter
        self.maxiter = maxiter
        self.additions = additions
        self.additionMap = additionMap
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        files = ["auc.stat", "color.map", "down.features", "include.samples",
                 "mut.circle", "mut.features", "mut.test.scores", "mut.train.scores",
                 "non.test.scores", "non.train.scores", "null.scores",
                 "paradigm_down.tab", "paradigm_up.tab", "pshift.train.tab",
                 "pshift.test.tab", "real.scores", "shift.circle", "sig.tab",
                 "up.features"]
        
        aucDiffs = {}
        f = open("iter_%s/auc.stat" % (self.iter-1), "r")
        try:
            lastAUC = float(f.readline().rstrip("\r\n").split("\t")[0])
        except ValueError:
            lastAUC = 0.0
        f.close()
        for addition in self.additionMap:
            if os.path.exists("iter_%s-%s/auc.stat" % (self.iter, addition)):
                f = open("iter_%s-%s/auc.stat" % (self.iter, addition), "r")
                try:
                    currentAUC = float(f.readline().rstrip("\r\n").split("\t")[0])
                except ValueError:
                    currentAUC = 0.0
                f.close()
                aucDiffs[addition] = currentAUC - lastAUC
        rankedAdditions = aucDiffs.keys()
        rankedAdditions.sort(lambda x, y: cmp(aucDiffs[y], aucDiffs[x]))
        if len(rankedAdditions) == 0:
            for file in files:
                if os.path.exists("iter_%s/%s" % (self.iter-1, file)):
                    system("cp iter_%s/%s %s" % (self.iter-1, file, file))
        elif aucDiffs[rankedAdditions[0]] <= 0.0:
            for file in files:
                if os.path.exists("iter_%s/%s" % (self.iter-1, file)):
                    system("cp iter_%s/%s %s" % (self.iter-1, file, file))        
        else:
            updatedAdditions = deepcopy(self.additions)
            updatedAdditions.append(rankedAdditions[0])
            system("ln -s iter_%s-%s iter_%s" % (self.iter, rankedAdditions[0], self.iter))
            if self.iter >= self.maxiter:
                for file in files:
                    if os.path.exists("iter_%s/%s" % (self.iter, file)):
                        system("cp iter_%s/%s %s" % (self.iter, file, file))
            elif len(set(self.additionMap.keys()) - set(updatedAdditions)) == 0:
                for file in files:
                    if os.path.exists("iter_%s/%s" % (self.iter, file)):
                        system("cp iter_%s/%s %s" % (self.iter, file, file))
            else:
                if self.fold == 0:
                    f = open("../../%s.log" % (self.mutatedFeature.gene), "a")
                    f.write("%s\t%s\t%s\n" % (self.additionMap[rankedAdditions[0]]))
                    f.close()
                self.setFollowOnTarget(branchLinkSearch(self.fold, self.mutatedFeature,
                                                          self.trainSamples, self.paradigmSetup,
                                                          self.uPathway, self.dPathway,
                                                          self.uBase, self.dBase,
                                                          self.extensions,
                                                          self.iter+1, self.maxiter,
                                                          updatedAdditions, self.additionMap,
                                                          self.directory))

class runPARADIGM(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, uPathway, 
                 dPathway, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## run paradigm (observed and nulls)
        system("echo Running PARADIGM inference ... >> progress.log")
        batches = batchCount(len(self.paradigmSetup.samples),
                             batchSize = self.paradigmSetup.size)
        if batches == 0:
            self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_" % (self.mutatedFeature.data), "%s_upstream.fa" % (self.mutatedFeature.gene)), self.directory))
            self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_" % (self.mutatedFeature.data), "%s_downstream.fa" % (self.mutatedFeature.gene)), self.directory))
            for null in range(1, self.paradigmSetup.nulls+1):
                self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_N%s_" % (self.mutatedFeature.data, null), "N%s_%s_upstream.fa" % (null, self.mutatedFeature.gene)), self.directory))
                self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_N%s_" % (self.mutatedFeature.data, null), "N%s_%s_downstream.fa" % (null, self.mutatedFeature.gene)), self.directory))
        elif self.paradigmSetup.public:
            system("mkdir outputFiles")
            for b in range(batches):
                self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_b%s_%s_" % (self.mutatedFeature.data, b, batches), "outputFiles/%s_upstream_b%s_%s.fa" % (self.mutatedFeature.gene, b, batches)), self.directory))
                self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_b%s_%s_" % (self.mutatedFeature.data, b, batches), "outputFiles/%s_downstream_b%s_%s.fa" % (self.mutatedFeature.gene, b, batches)), self.directory))
                for null in range(1, self.paradigmSetup.nulls+1):
                    self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/up_N%s_b%s_%s_" % (self.mutatedFeature.data, null, b, batches), "outputFiles/N%s_%s_upstream_b%s_%s.fa" % (null, self.mutatedFeature.gene, b, batches)), self.directory))
                    self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s" % (paradigmExec, "%s/down_N%s_b%s_%s_" % (self.mutatedFeature.data, null, b, batches), "outputFiles/N%s_%s_downstream_b%s_%s.fa" % (null, self.mutatedFeature.gene, b, batches)), self.directory))
        else:
            system("mkdir outputFiles")
            for b in range(batches):
                self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/up_" % (self.mutatedFeature.data), "outputFiles/%s_upstream_b%s_%s.fa" % (self.mutatedFeature.gene, b, batches), b, batches), self.directory))
                self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/down_" % (self.mutatedFeature.data), "outputFiles/%s_downstream_b%s_%s.fa" % (self.mutatedFeature.gene, b, batches), b, batches), self.directory))
                for null in range(1, self.paradigmSetup.nulls+1):
                    self.addChildTarget(jtCmd("%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/up_N%s_" % (self.mutatedFeature.data, null), "outputFiles/N%s_%s_upstream_b%s_%s.fa" % (null, self.mutatedFeature.gene, b, batches), b, batches), self.directory))
                    self.addChildTarget(jtCmd("%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s" % (paradigmExec, "%s/down_N%s_" % (self.mutatedFeature.data, null), "outputFiles/N%s_%s_downstream_b%s_%s.fa" % (null, self.mutatedFeature.gene, b, batches), b, batches), self.directory))
        self.setFollowOnTarget(collectPARADIGM(self.fold, self.mutatedFeature, 
                                          self.trainSamples, self.paradigmSetup,
                                          self.uPathway, self.dPathway, self.directory))

class collectPARADIGM(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, uPathway, 
                 dPathway, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        batches = batchCount(len(self.paradigmSetup.samples),
                             batchSize = self.paradigmSetup.size)
        for b in range(batches):
            system("cat outputFiles/%s_upstream_b%s_%s.fa >> %s_upstream.fa" % (self.mutatedFeature.gene, b, batches, self.mutatedFeature.gene))
            system("cat outputFiles/%s_downstream_b%s_%s.fa >> %s_downstream.fa" % (self.mutatedFeature.gene, b, batches, self.mutatedFeature.gene))
            for null in range(1, self.paradigmSetup.nulls+1):
                system("cat outputFiles/N%s_%s_upstream_b%s_%s.fa >> N%s_%s_upstream.fa" % (null, self.mutatedFeature.gene, b, batches, null, self.mutatedFeature.gene))
                system("cat outputFiles/N%s_%s_downstream_b%s_%s.fa >> N%s_%s_downstream.fa" % (null, self.mutatedFeature.gene, b, batches, null, self.mutatedFeature.gene))
        if os.path.exists("outputFiles"):
            system("rm -rf outputFiles")
        
        self.setFollowOnTarget(pshiftCV(self.fold, self.mutatedFeature, self.trainSamples, self.paradigmSetup, self.uPathway, self.dPathway, self.directory, nulls = self.paradigmSetup.nulls))

class pshiftCV(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, uPathway, 
                 dPathway, directory, nulls = 10, msepMethod = "tt", alpha = 0.05):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.uPathway = uPathway
        self.dPathway = dPathway
        self.nulls = nulls
        self.msepMethod = msepMethod
        self.alpha = alpha
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## define sample groups
        trainGroup = self.trainSamples
        testGroup = list((set(self.paradigmSetup.samples) - set(self.trainSamples)) & 
             (set(self.mutatedFeature.positive) | set(self.mutatedFeature.negative)))
        mutGroup_tr = list(set(self.mutatedFeature.positive) & set(trainGroup))
        mutGroup_te = list(set(self.mutatedFeature.positive) & set(testGroup))
        nonGroup_tr = list(set(self.mutatedFeature.negative) & set(trainGroup))
        nonGroup_te = list(set(self.mutatedFeature.negative) & set(testGroup))
        
        ## read paradigm output
        upScore = rPARADIGM("%s_upstream.fa" % (self.mutatedFeature.gene), 
                            useRows = self.uPathway.nodes.keys())[1]
        downScore = rPARADIGM("%s_downstream.fa" % (self.mutatedFeature.gene),
                            useRows = self.dPathway.nodes.keys())[1]
        wCRSData("paradigm_up.tab", upScore)
        wCRSData("paradigm_down.tab", downScore)
        n_upScore = {}
        n_downScore = {}
        for null in range(1, self.nulls+1):
            if os.path.exists("N%s_%s_upstream.fa" % (null, self.mutatedFeature.gene)):
                n_upScore[null] = rPARADIGM("N%s_%s_upstream.fa" %
                                            (null, self.mutatedFeature.gene),
                                            useRows = [self.mutatedFeature.gene])[1]
            if os.path.exists("N%s_%s_downstream.fa" % (null, self.mutatedFeature.gene)):
                n_downScore[null] = rPARADIGM("N%s_%s_downstream.fa" %
                                              (null, self.mutatedFeature.gene),
                                              useRows = [self.mutatedFeature.gene])[1]
    
        ## score raw shifts
        rawShift = {}
        for sample in self.paradigmSetup.samples:
            if (sample in downScore) and (sample in upScore):
                rawShift[sample] = (downScore[sample][self.mutatedFeature.gene] -
                                     upScore[sample][self.mutatedFeature.gene])
        for sample in self.paradigmSetup.samples:
            if (sample in downScore) and (sample in upScore):
                for null in range(1, self.nulls+1):
                    rawShift["null%s_%s" % (null, sample)] = (
                                     n_downScore[null][sample][self.mutatedFeature.gene] -
                                     n_upScore[null][sample][self.mutatedFeature.gene])
        
        ## store shifts for mutants and non-mutants normalized by trained non-mutants
        (nonMean, nonStd) = mean_std([rawShift[sample] for sample in nonGroup_tr])
        labelMap = {}
        normShift_tr = {}
        normShift_te = {}
        o = open("pshift.train.tab", "w")
        o.write("> %s\tP-Shifts:Table\n" % (self.mutatedFeature.gene))
        o.write("# sample\tclass\tP-Shift\n")
        f = open("mut.train.scores", "w")
        for sample in mutGroup_tr:
            labelMap[sample] = 1
            normShift_tr[sample] = (rawShift[sample]-nonMean)
            o.write("%s\t+\t%s\n" % (sample, normShift_tr[sample]))
            f.write("%s\t%s\n" % (sample, normShift_tr[sample]))
        f.close()
        f = open("non.train.scores", "w")
        for sample in nonGroup_tr:
            labelMap[sample] = 0
            normShift_tr[sample] = (rawShift[sample]-nonMean)
            o.write("%s\t-\t%s\n" % (sample, normShift_tr[sample]))
            f.write("%s\t%s\n" % (sample, normShift_tr[sample]))
        f.close()
        o.close()
        o = open("pshift.test.tab", "w")
        o.write("> %s\tP-Shifts:Table\n" % (self.mutatedFeature.gene))
        o.write("# sample\tclass\tP-Shift\n")
        f = open("mut.test.scores", "w")
        for sample in mutGroup_te:
            labelMap[sample] = 1
            normShift_te[sample] = (rawShift[sample]-nonMean)
            o.write("%s\t+\t%s\n" % (sample, normShift_te[sample]))
            f.write("%s\t%s\n" % (sample, normShift_te[sample]))
        f.close()
        f = open("non.test.scores", "w")
        for sample in nonGroup_te:
            labelMap[sample] = 0
            normShift_te[sample] = (rawShift[sample]-nonMean)
            o.write("%s\t-\t%s\n" % (sample, normShift_te[sample]))
            f.write("%s\t%s\n" % (sample, normShift_te[sample]))
        f.close()
        o.close()
    
        ## compute auc from stats
        sorted_tr = normShift_tr.keys()
        sorted_tr.sort(lambda x, y: cmp(normShift_tr[y],normShift_tr[x]))
        auc_tr = computeAUC(sorted_tr, labelMap, normShift_tr)[0]
        if self.fold != 0:
            sorted_te = normShift_te.keys()
            sorted_te.sort(lambda x, y: cmp(normShift_te[y],normShift_te[x]))
            auc_te = computeAUC(sorted_te, labelMap, normShift_te)[0]
        else:
            auc_te = "---"
        f = open("auc.stat", "w")
        f.write("%s\t%s\n" % (auc_tr, auc_te))
        f.close()
        
        ## output files for circle maps
        if self.fold == 0:
            f = open("up.features", "w")
            for feature in (set(upScore[upScore.keys()[0]].keys()) -
                                    set([self.mutatedFeature.gene])):
                f.write("%s\n" % (feature))
            f.close()
            f = open("down.features", "w")
            for feature in (set(downScore[downScore.keys()[0]].keys()) - 
                                        set([self.mutatedFeature.gene])):
                f.write("%s\n" % (feature))
            f.close()
            f = open("mut.features", "w")
            f.write("%s\n" % (self.mutatedFeature.gene))
            f.close()
            f = open("include.samples", "w")
            for sample in (self.paradigmSetup.samples):
                f.write("%s\n" % (sample))
            f.close()
            f = open("mut.circle", "w")
            f.write("id")
            for sample in (self.paradigmSetup.samples):
                if sample not in rawShift:
                    continue
                f.write("\t%s" % (sample))
            f.write("\n")
            f.write("*")
            for sample in (self.paradigmSetup.samples):
                if sample not in rawShift:
                    continue
                if sample in self.mutatedFeature.positive:
                    f.write("\t1")
                elif sample in self.mutatedFeature.negative:
                    f.write("\t0")
                else:
                    f.write("\t0.5")
            f.write("\n")
            f.close()
            f = open("shift.circle", "w")
            f.write("id")
            for sample in (self.paradigmSetup.samples):
                if sample not in rawShift:
                    continue
                f.write("\t%s" % (sample))
            f.write("\n")
            f.write("*")
            for sample in (self.paradigmSetup.samples):
                if sample not in rawShift:
                    continue
                f.write("\t%s" % (rawShift[sample]))
            f.write("\n")
            f.close()
            f = open("color.map", "w")
            f.write("> 1\n0\t255.255.255\n1\t0.0.0\n0.5\t150.150.150\n")
            f.close()

        ## compute mutant separation
        if self.fold == 0:
            msepScore = {}
            msepScore["real"] = mutSeparation(nonGroup_tr, mutGroup_tr, rawShift, 
                                              method = self.msepMethod)
            for null in range(1, self.nulls+1):
                msepScore["null%s" % (null)] = mutSeparation(["null%s_%s" % (null, sample) 
                                                             for sample in nonGroup_tr],
                                                             ["null%s_%s" % (null, sample)
                                                             for sample in mutGroup_tr],
                                                             rawShift, 
                                                             method = self.msepMethod)
            f = open("real.scores", "w")
            f.write("%s\n" % (msepScore["real"]))
            f.close()
            f = open("null.scores", "w")
            for null in range(1, self.nulls+1):
                f.write("%s\n" % (msepScore["null%s" % (null)]))
            f.close()
    
            ## compute background significance
            realScores = [msepScore["real"]]
            nullScores = [msepScore[null] for null in list(set(msepScore.keys()) - set(["real"]))]
            (nullMean, nullStd) = mean_std(nullScores)
            significanceScore = (realScores[0]-nullMean)/(nullStd+self.alpha)
        
            ## output sig.stats
            f = open("sig.tab", "w")
            f.write("# gene\tnon_mut\tmut\tmsep\tsignificance\n")
            f.write("%s\t%s\t%s\t%s\t%s\n" % (self.mutatedFeature.gene, len(nonGroup_tr),
                                              len(mutGroup_tr), 
                                              msepScore["real"], significanceScore))
            f.close()
        else:
            f = open("sig.tab", "w")
            f.write("# gene\tnon_mut\tmut\tmsep\tsignificance\n")
            f.write("%s\t%s\t%s\t%s\t%s\n" % (self.mutatedFeature.gene,
                                                  len(nonGroup_tr), len(mutGroup_tr),
                                                  "---", "---"))
            f.close()

class evaluateParams(Target):
    def __init__(self, fold, mutatedFeature, trainSamples, paradigmSetup, gPathway, params, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.mutatedFeature = mutatedFeature
        self.trainSamples = trainSamples
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.params = params
        self.directory = directory
    def run(self):
        if self.fold == 0:
            shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedFeature.gene)
            os.chdir(shiftDir)
            system("echo Identifying Best Param by Training Validation ... >> progress.log")
        else:
            shiftDir = "%s/analysis/%s/fold%s" % (self.directory, self.mutatedFeature.gene, self.fold)
            os.chdir(shiftDir)
        
        ## report AUCs for model with best training validation
        topAUC_tr = 0
        topAUC_te = 0
        topAUC_params = None
        for distance in self.params.distance:
            for thresh in self.params.threshold:
                for inc in self.params.penalty:
                    for method in self.params.selection:
                        f = open("param_%s_%s_%s_%s/auc.stat" % (distance, thresh, inc, method), "r")
                        line = f.readline()
                        f.close()
                        (auc_tr, auc_te) = line.rstrip().split("\t")
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
                            topAUC_params = [str(distance), str(thresh), str(inc), str(method)]
        if self.fold != 0:
            f = open("auc.stat", "w")
            if topAUC_tr == 0:
                f.write("NA\tNA\tNA\n")
            else:
                f.write("%s\t%s\t%s\n" % (topAUC_tr, topAUC_te, ",".join(topAUC_params)))
            f.close()
        else:
            if topAUC_params is not None:
                self.setFollowOnTarget(mutationSHIFT(topAUC_params, self.mutatedFeature, self.paradigmSetup, self.gPathway, self.directory, clean = self.params.clean))
            
class mutationSHIFT(Target):
    def __init__(self, topAUC_params, mutatedFeature, paradigmSetup, gPathway, directory, clean = True):
        Target.__init__(self, time=10000)
        self.topAUC_params = topAUC_params
        self.mutatedFeature = mutatedFeature
        self.paradigmSetup = paradigmSetup
        self.gPathway = gPathway
        self.directory = directory
        self.clean = clean
    def run(self):
        shiftDir = "%s/analysis/%s" % (self.directory, self.mutatedFeature.gene)
        os.chdir(shiftDir)
        
        ## copy tables
        system("echo Output Reports ... >> progress.log")
        system("cp param_%s/sig.tab sig.tab" % ("_".join(self.topAUC_params)))
        system("cp param_%s/pshift.train.tab pshift.tab" % ("_".join(self.topAUC_params)))
        
        ## output msep and background plots
        system("%s/msep.R %s param_%s/mut.train.scores param_%s/non.train.scores" % (self.directory, 
                                                 self.mutatedFeature.gene, "_".join(self.topAUC_params), 
                                                 "_".join(self.topAUC_params)))
        system("%s/background.R %s param_%s/real.scores param_%s/null.scores" % (self.directory, 
                                                 self.mutatedFeature.gene, "_".join(self.topAUC_params), 
                                                 "_".join(self.topAUC_params)))
        
        ## prepare web report
        system("mkdir img")
        
        ## draw circles
        system("cat data/up_%s | transpose.pl > param_%s/exp.tab" % (self.paradigmSetup.exp.split("/")[-1], "_".join(self.topAUC_params)))
        os.chdir("param_%s" % ("_".join(self.topAUC_params)))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f mut.features ../img/ mut.circle exp.tab paradigm_up.tab paradigm_down.tab shift.circle" % (circleExec, self.mutatedFeature.gene))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f up.features ../img/ mut.circle exp.tab paradigm_up.tab" % (circleExec, self.mutatedFeature.gene))
        system("%s -o \"%s;mut.circle,shift.circle\" -s include.samples -f down.features ../img/ mut.circle exp.tab paradigm_down.tab" % (circleExec, self.mutatedFeature.gene))
        os.chdir("..")
        
        ## output sif
        (uNodes, uInteractions) = rPathway("param_%s/upstream_pathway.tab" %
                                                    ("_".join(self.topAUC_params)))
        (dNodes, dInteractions) = rPathway("param_%s/downstream_pathway.tab" % 
                                                    ("_".join(self.topAUC_params)))
        cPathway = combinePathways(Pathway(uNodes, uInteractions),
                                   Pathway(dNodes, dInteractions))
        wSIF("pshift_%s.sif" % (self.mutatedFeature.gene), cPathway.interactions)
        
        ## node attributes
        scoreMap = {}
        scoreMap["pshift_%s" % (self.mutatedFeature.gene)] = {}
        for node in self.gPathway.nodes.keys():
            if node == self.mutatedFeature.gene:
                scoreMap["pshift_%s" % (self.mutatedFeature.gene)][node] = 10
            else:
                scoreMap["pshift_%s" % (self.mutatedFeature.gene)][node] = 7
        wNodeAttributes(self.gPathway.nodes, scoreMap = scoreMap, directory = "./")
        
        ## clean
        if self.clean:
            system("rm -rf data")

class pshiftReport(Target):
    def __init__(self, includeFeatures, reportDir, directory):
        Target.__init__(self, time=10000)
        self.includeFeatures = includeFeatures
        self.reportDir = reportDir
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## cytoscape-web
        for gene in self.includeFeatures:
            if os.path.exists("analysis/%s/sig.tab" % (gene)):
                tableFiles = []
                tableFiles.append("analysis/%s/sig.tab" % (gene))
                tableFiles.append("msepPlot:analysis/%s/%s.msep.pdf" % (gene, gene))
                tableFiles.append("backgroundPlot:analysis/%s/%s.background.pdf" %
                                                                             (gene, gene))
                tableFiles.append("analysis/%s/avgAUC.tab" % (gene))
                tableFiles.append("analysis/%s/pshift.tab" % (gene))
                system("pathmark-report.py -t %s analysis/%s %s" %
                                             (",".join(tableFiles), gene, self.reportDir))
                system("cp analysis/%s/pshift* %s" % (gene, self.reportDir))

def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        print "WARNING: .jobTree directory detected, remove it first to start a fresh run"
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] paradigmDir mutFile ")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile", help="Add as child of jobFile rather than new jobTree")
    parser.add_option("-c", "--config", dest="configFile", default=None)
    parser.add_option("-s", "--samples", dest="includeSamples", default=None)
    parser.add_option("-f", "--features", dest="includeFeatures", default=None)
    parser.add_option("-p", "--public", action="store_true", dest="paradigmPublic",
                      default=False)
    parser.add_option("-b", "--batchSize", dest="batchSize", default=50)
    parser.add_option("-n", "--nulls", dest="nulls", default=30)
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
    
    ## get params
    params = Parameters()
    if options.configFile:
        params.set(options.configFile)
    
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
    
    ## get paradigm files
    paradigmSetup = ParadigmSetup(paradigmDir, options.includeSamples,
                                  public = options.paradigmPublic,
                                  size = int(options.batchSize),
                                  nulls = int(options.nulls))
    
    ## get pathway
    (gNodes, gInteractions) = rPathway(paradigmSetup.pathway)
    gPathway = Pathway(gNodes, gInteractions)
    
    ## get mutation annotations
    mutatedList = []
    mutatedMap = {}
    f = open(mutFile, "r")
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split("\t")
        if len(pline) == 3:
            mutatedFeature = Mutated(pline[0], paradigmSetup.samples, pline[2].split(","),
                                     negative = None, directory = os.getcwd())
        elif len(pline) == 4:
            mutatedFeature = Mutated(pline[0], paradigmSetup.samples, pline[2].split(","),
                                     negative = pline[3].split(","),
                                     directory = os.getcwd())
        if mutatedFeature.gene in gPathway.nodes:
            if len(mutatedFeature.positive) >= mutationThreshold:
                mutatedMap[mutatedFeature.gene] = mutatedFeature
                mutatedList.append(mutatedFeature.gene)
    f.close()
    if options.includeFeatures:
        mutatedList = []
        for gene in rList(options.includeFeatures):
            if gene in mutatedMap:
                mutatedList.append(gene)
    submitList = []
    for gene in mutatedList:
        submitList.append(deepcopy(mutatedMap[gene]))
        if len(submitList) >= maxFeatures:
            break
    
    ## run
    logger.info("options: " + str(options))
    logger.info("starting make")
    writeScripts()
    if runParallel:
        s = Stack(branchGenes(submitList, paradigmSetup, gPathway, params, foldMap,
                              os.getcwd()))
    else:
        s = Stack(queueGenes(submitList, paradigmSetup, gPathway, params, foldMap,
                             os.getcwd()))
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
