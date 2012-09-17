import math, os, sys, random, re
from copy import deepcopy
import liblinear
import liblinearutil

import mData, mCalculate

def log(msg, die = False):
    """logger function"""
    sys.stderr.write(msg)
    if die:
        sys.exit(1)

class Pathway:
    def __init__(self, nodes, interactions, pid = None):
        self.nodes = nodes
        self.interactions = interactions
        self.pid = pid
    def removeNode(self, node):
        (self.nodes, self.interactions) = removeNode(node, self.nodes, self.interactions)
    def selfTest(self):
        ## check for unannotated features
        for source in self.interactions.keys():
            if source not in self.nodes:
                if source.endswith("(complex)"):
                    log("WARNING: %s not defined, default to complex\n" % (source))
                    self.nodes[source] = "complex"
                elif source.endswith("(family)"):
                    log("WARNING: %s not defined, default to family\n" % (source))
                    self.nodes[source] = "family"
                elif source.endswith("(abstract)"):
                    log("WARNING: %s not defined, default to abstract\n" % (source))
                    self.nodes[source] = "abstract"
                else:
                    log("WARNING: %s not defined, default to protein\n" % (source))
                    self.nodes[source] = "protein"
            for target in self.interactions[source].keys():
                if target not in self.nodes:
                    if target.endswith("(complex)"):
                        log("WARNING: %s not defined, default to complex\n" % (target))
                        self.nodes[target] = "complex"
                    elif target.endswith("(family)"):
                        log("WARNING: %s not defined, default to family\n" % (target))
                        self.nodes[target] = "family"
                    elif target.endswith("(abstract)"):
                        log("WARNING: %s not defined, default to abstract\n" % (target))
                        self.nodes[target] = "abstract"
                    else:
                        log("WARNING: %s not defined, default to protein\n" % (target))
                        self.nodes[target] = "protein"
        ## check for invalid links
        legalLinks = ["-t>", "-t|", "-a>", "-a|", "-ap>", "-ap|", "component>", "member>"]
        for source in self.interactions.keys():
            for target in self.interactions[source].keys():
                for link in re.split(";", self.interactions[source][target]):
                    if link not in legalLinks:
                        log("WARNING: illegal link %s\t%s\t%s found\n" % (source, target, link))
        ## report number of separate components
        components = sortConnected(self)
        log("found %s subpathway components\n" % (len(components)))
        (self.nodes, self.interactions) = constructInteractions(components[0], self.nodes, self.interactions)
    def wAliasMap(self, hugof, outf):
        hugoList = mData.rList(hugof)
        if self.pid != None:
            prefix = "%s_" % (self.pid)
        else:
            prefix = ""
        f = open(outf, "w")
        for node in (set(hugoList) & set(self.nodes.keys())):
            f.write("%s\t%s\n" % (prefix+node, node))
        f.close()

def rPathway(inf, reverse = False, retProteins = False, delim = "\t"):
    """read UCSC pathway format"""
    proteins = set()
    readPathway = Pathway(dict(), dict())
    f = open(inf, "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if len(pline) == 2:
            readPathway.nodes[pline[1]] = pline[0]
            if pline[0] == "protein":
                proteins.update([pline[1]])
        elif len(pline) == 3:
            if reverse:
                if pline[1] not in readPathway.interactions:
                    readPathway.interactions[pline[1]] = dict()
                if pline[0] not in readPathway.interactions[pline[1]]:
                    readPathway.interactions[pline[1]][pline[0]] = pline[2]
                else:
                    readPathway.interactions[pline[1]][pline[0]] += ";"+pline[2]
            else:
                if pline[0] not in readPathway.interactions:
                    readPathway.interactions[pline[0]] = dict()
                if pline[1] not in readPathway.interactions[pline[0]]:
                    readPathway.interactions[pline[0]][pline[1]] = pline[2]
                else:
                    readPathway.interactions[pline[0]][pline[1]] += ";"+pline[2]
        else:
            print >> sys.stderr, "ERROR: line length not 2 or 3: \"%s\"" % (line)
            sys.exit(1)
    f.close()
    if retProteins:
        return(readPathway.nodes, readPathway.interactions, proteins)
    else:
        return(readPathway.nodes, readPathway.interactions)

def wPathway(outf, outNodes, outInteractions, useNodes = None):
    """write UCSC pathway format"""
    f = open(outf, "w")
    if useNodes == None:
        useNodes = outNodes.keys()
    for i in useNodes:
        if i not in outNodes:
            continue
        f.write("%s\t%s\n" % (outNodes[i], i))
    for i in useNodes:
        if i not in outInteractions:
            continue
        for j in outInteractions[i].keys():
            if j not in useNodes:
                continue
            for k in re.split(";", outInteractions[i][j]):
                f.write("%s\t%s\t%s\n" % (i, j, k))
    f.close()

def wSIF(writeFile, writeInteractions, useNodes = None):
    """write .sif"""
    f = open(writeFile, "w")
    if useNodes == None:
        for i in writeInteractions.keys():
            for j in writeInteractions[i].keys():
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    else:
        for i in useNodes:
            if i not in writeInteractions:
                continue
            for j in writeInteractions[i].keys():
                if j not in useNodes:
                    continue
                for k in re.split(";", writeInteractions[i][j]):
                    f.write("%s\t%s\t%s\n" % (i, k, j))
    f.close()

def reverseInteractions(interactions):
    """reverse interaction mapping"""
    rinteractions = {}
    for source in interactions.keys():
        for target in interactions[source].keys():
            if target not in rinteractions:
                rinteractions[target] = {}
            rinteractions[target][source] = interactions[source][target]
    return(rinteractions)

def getComponentMap(pNodes, pInteractions):
    """create the dictionary componentMap from interaction map"""
    rpInteractions = reverseInteractions(pInteractions)
    componentMap = dict()
    for i in pNodes.keys():
        if pNodes[i] != "complex":
            continue
        componentMap[i] = []
        if i not in rpInteractions:
            continue
        for j in rpInteractions[i]:
            if rpInteractions[i][j] == "component>":
                componentMap[i].append(j)
    return(componentMap)

def constructInteractions(nodeList, refNodes, refInteractions):
    """select concepts from list and construct Pathway"""
    outPathway = Pathway({}, {})
    for i in nodeList:
        outPathway.nodes[i] = refNodes[i]
        if i in refInteractions:
            for j in refInteractions[i].keys():
                if j in nodeList:
                    if i not in outPathway.interactions:
                        outPathway.interactions[i] = dict()
                    outPathway.interactions[i][j] = refInteractions[i][j]
    return(outPathway.nodes, outPathway.interactions)

def combinePathways(currentPathway, appendPathway, exclude = []):
    """combine Pathways"""
    for source in appendPathway.interactions.keys():
        if source in exclude:
            continue
        if source not in currentPathway.nodes:
            currentPathway.nodes[source] = appendPathway.nodes[source]
        for target in appendPathway.interactions[source].keys():
            if target in exclude:
                continue
            if target not in currentPathway.nodes:
                currentPathway.nodes[target] = appendPathway.nodes[target]
            if source not in currentPathway.interactions:
                currentPathway.interactions[source] = dict()
            if target not in currentPathway.interactions[source]:
                currentPathway.interactions[source][target] = appendPathway.interactions[source][target]
    return(currentPathway)

def sortConnected(inPathway):
    rinteractions = reverseInteractions(inPathway.interactions)
    index = 1
    mapNets = dict()
    sortedNets = []
    seenNodes = set()
    for i in inPathway.nodes.keys():
        if i in seenNodes:
            continue
        borderNodes = [i]
        currentNet = [i]
        while len(borderNodes) > 0:
            if borderNodes[0] in rinteractions:
                for j in rinteractions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            if borderNodes[0] in inPathway.interactions:
                for j in inPathway.interactions[borderNodes[0]].keys():
                    if j not in seenNodes:
                        seenNodes.update([j])
                        borderNodes.append(j)
                        currentNet.append(j)
            borderNodes.pop(0)
        if ("__DISCONNECTED__" not in currentNet):
            mapNets[index] = deepcopy(currentNet)
            index += 1
    indexList = mapNets.keys()
    netScore = dict()
    for i in indexList:
        netScore[i] = len(mapNets[i])
    indexList.sort(lambda x, y: cmp(netScore[y], netScore[x]))
    for i in indexList:
        sortedNets.append(mapNets[i])
    return(sortedNets)

def flattenPathway(inPathway):
    """expands complexes into their respective gene components"""
    allowedNodes = ["abstract", "family", "miRNA", "protein", "rna"]
    outPathway = Pathway({}, {})
    ## read and search componentMap for protein components
    componentMap = getComponentMap(inPathway.nodes, inPathway.interactions)
    for entity in componentMap.keys():
        seenNodes = set()
        elements = []
        expand = deepcopy(componentMap[entity])
        while len(expand) > 0:
            if expand[0] in seenNodes:
                expand.pop(0)
                continue
            seenNodes.update([expand[0]])
            if inPathway.nodes[expand[0]] == "protein":
                elements.append(expand[0])
            elif expand[0] in componentMap:
                expand += deepcopy(componentMap[expand[0]])
            expand.pop(0)
        componentMap[entity] = elements
    ## iterate over all interactions
    for source in inPathway.interactions.keys(): 
        for target in inPathway.interactions[source].keys():
            ## update interactions map
            if inPathway.nodes[source] in allowedNodes:
                if inPathway.nodes[target] in allowedNodes:
                    if inPathway.interactions[source][target] not in ["component>", "member>"]:
                        if source not in outPathway.nodes:
                            outPathway.nodes[source] = inPathway.nodes[source]
                        if target not in outPathway.nodes:
                            outPathway.nodes[target] = inPathway.nodes[target]
                        if source not in outPathway.interactions:
                            outPathway.interactions[source] = {}
                        outPathway.interactions[source][target] = inPathway.interactions[source][target]
                elif target in componentMap:
                    for element in componentMap[target]:
                        if source != element:
                            if source not in outPathway.nodes:
                                outPathway.nodes[source] = inPathway.nodes[source]
                            if element not in outPathway.nodes:
                                outPathway.nodes[element] = inPathway.nodes[element]
                            if source not in outPathway.interactions:
                                outPathway.interactions[source] = {}
                            if inPathway.interactions[source][target] == "component>":
                                outPathway.interactions[source][element] = "-a>"
                            else:
                                outPathway.interactions[source][element] = inPathway.interactions[source][target]
            elif source in componentMap:
                if inPathway.nodes[target] in allowedNodes:
                    for element in componentMap[source]:
                        if element not in outPathway.nodes:
                            outPathway.nodes[element] = inPathway.nodes[element]
                        if target not in outPathway.nodes:
                            outPathway.nodes[target] = inPathway.nodes[target]
                        if element not in outPathway.interactions:
                            outPathway.interactions[element] = {}
                        outPathway.interactions[element][target] = inPathway.interactions[source][target]
                elif target in componentMap:
                    continue
    return(outPathway)

def getNeighbors(node, distance, interactions):
    """returns upstream and downstream neighbors of distance"""
    rinteractions = reverseInteractions(interactions)
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for d in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in interactions:
                for i in interactions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
            if currNode in rinteractions:
                for i in rinteractions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = []
    return(seenNodes)

def getDownstream(node, distance, interactions):
    """returns downstream neighbors of distance"""
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for dist in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in interactions:
                for i in interactions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = list()
    return(list(seenNodes))

def getUpstream(node, distance, interactions):
    """returns downstream neighbors of distance"""
    rinteractions = reverseInteractions(interactions)
    seenNodes = set([node])
    borderNodes = [node]
    frontierNodes = []
    for dist in range(distance):
        while len(borderNodes) > 0:
            currNode = borderNodes.pop()
            if currNode in rinteractions:
                for i in rinteractions[currNode].keys():
                    if i not in seenNodes:
                        seenNodes.update([i])
                        frontierNodes.append(i)
        borderNodes = deepcopy(frontierNodes)
        frontierNodes = []
    return(list(seenNodes))

def shortestPath(source, target, interactions):
    shortestPaths = []
    allPaths = [[source]]
    while len(shortestPaths) == 0:
        nextPaths = []
        while len(allPaths) > 0:
            currentPath = allPaths.pop()
            if currentPath[-1] not in interactions:
                continue
            for intermediate in interactions[currentPath[-1]]:
                if intermediate in currentPath:
                    continue
                if intermediate == target:
                    shortestPaths.append(currentPath + [intermediate])
                    nextPaths.append(currentPath + [intermediate])
                else:
                    nextPaths.append(currentPath + [intermediate])
        if len(nextPaths) == 0:
            log("ERROR: could not find path between %s and %s\n" % (source, target), die = True)
        allPaths = deepcopy(nextPaths)
    return (shortestPaths)

def wNodeAttributes(pNodes, scoreMap = None, directory = "."):
    """write cytoscape node attribute files"""
    ## TYPE.NA
    typef = open("%s/TYPE.NA" % (directory), "w")
    typef.write("TYPE (class=java.lang.String)\n")
    for node in pNodes.keys():
        typef.write("%s = %s\n" % (node, pNodes[node]))
    typef.close()
    ## LABEL.NA
    labelf = open("%s/LABEL.NA" % (directory), "w")
    labelf.write("LABEL (class=java.lang.String)\n")
    for node in pNodes.keys():
        if pNodes[node] == "protein":
            labelf.write("%s = %s\n" % (node, node))
        else:
            labelf.write("%s = %s\n" % (node, ""))
    labelf.close()
    ## SCORE.NA
    if scoreMap != None:
        for element in scoreMap.keys():
            scoref = open("%s/SCORE.NA" % (directory), "w")
            scoref.write("SCORE (class=java.lang.Double)\n")
            for node in scoreMap[element].keys():
                scoref.write("%s = %s\n" % (node, scoreMap[element][node]))
            scoref.close()
   
def removeNode(node, pNodes, pInteractions):
    """remove a node and its interactions from a Pathway"""
    rpInteractions = reverseInteractions(pInteractions)
    del pNodes[node]
    if node in pInteractions:
        for element in pInteractions[node].keys():
            del pInteractions[node][element]
            if len(pInteractions[node].keys()) == 0:
                del pInteractions[node]
            del rpInteractions[element][node]
            if len(rpInteractions[element].keys()) == 0:
                del rpInteractions[element]
    if node in rpInteractions:
        for element in rpInteractions[node].keys():
            del pInteractions[element][node]
            if len(pInteractions[element].keys()) == 0:
                del pInteractions[element]
            del rpInteractions[node][element]
            if len(rpInteractions[node].keys()) == 0:
                del rpInteractions[node]
    return(pNodes, pInteractions)

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

def mean(inList, null = "NA"):
    """Calculates mean"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = null
    else:
        mean = sum(fList)/len(fList)
    return (mean)

def ttest(values0, values1, alpha = 0.05):
    (mean0, std0) = mCalculate.mean_std(values0)
    (mean1, std1) = mCalculate.mean_std(values1)
    try:
        tval = (mean1-mean0)/(math.sqrt((std1**2)/len(values1)+(std0**2)/len(values0))+alpha)
    except ZeroDivisionError:
        tval = "NA"
    except TypeError:
        tval = "NA"
    return(tval)

def kstest(values0, values1):
    (mean0, std0) = mCalculate.mean_std(values0)
    (mean1, std1) = mCalculate.mean_std(values1)
    values0.sort()
    values1.sort()
    ksval = 0.0
    if (len(values0) == 0) | (len(values1) == 0):
        return("NA")
    for i, y in enumerate(values1):
        chigh = sum([x <= y for x in values0])/float(len(values0)) - (i)/float(len(values1))
        clow = (i+1)/float(len(values1)) - sum([x < y for x in values0])/float(len(values0))
        if chigh >= clow:
            val = chigh
        else:
            val = -clow
        if abs(val) > abs(ksval):
            ksval = val
            point = y
    ks = open("ks.stat", "w")
    ks.write("ks.test\t%s\t%s\n" % (point, ksval))
    ks.close()
    return(ksval)

def mutSeparation(group0, group1, data, method = "tt"):
    """computes the signal score from shiftScore"""
    values0 = []
    values1 = []
    for i in (group0 + group1):
        if data[i] == "NA":
            continue
        if i in group0:
            values0.append(data[i])
        elif i in group1:
            values1.append(data[i])
    if method == "tt":
        val = ttest(values0, values1)
    elif method == "ks":
        val = kstest(values0, values1)
    else:
        val = "NA"
    return(val)

def absScores(scoreMap):
    absMap = {}
    for feature in scoreMap.keys():
        try:
            fval = float(scoreMap[feature])
            if fval != fval:
                raise ValueError
            absMap[feature] = abs(fval)
        except ValueError:
            absMap[feature] = 0
    return(absMap)

def rankScores(scoreMap):
    rankMap = {}
    features = []
    for feature in scoreMap.keys():
        try:
            fval = float(scoreMap[feature])
            if fval != fval:
                raise ValueError
            features.append(feature)
        except ValueError:
            rankMap[feature] = 0
    features.sort(lambda x, y: cmp(abs(scoreMap[x]),abs(scoreMap[y])))
    for i, feature in enumerate(features):
        rankMap[feature] = float(i+1)/len(features)
    return(rankMap)

def scoreSVM(matData, posSamples, negSamples):
    svmLabels = []
    svmData = []
    featureMap = {}
    for sample in (posSamples + negSamples):
        if sample in posSamples:
            svmLabels.append(1)
        else:
            svmLabels.append(-1)
        svmData.append({})
        for i, feature in enumerate(matData.keys()):
            if i+1 not in featureMap:
                featureMap[i+1] = feature
            try:
                svmData[-1][i+1] = float(matData[feature][sample])
            except ValueError:
                svmData[-1][i+1] = 0.0
    prob = liblinear.problem(svmLabels, svmData)
    param = liblinear.parameter('-s 3 -c 5 -q')
    liblinearutil.save_model('model_file', liblinearutil.train(prob, param))
    weights = mData.rList("model_file")[6:]
    scoreMap = {}
    for feature in featureMap.keys():
        scoreMap[featureMap[feature]] = float(weights[feature-1])
    return (scoreMap)

def scoreTT(matData, posSamples, negSamples):
    scoreMap = {}
    for feature in matData.keys():
        scoreMap[feature] = ttest([matData[feature][i] for i in negSamples], [matData[feature][i] for i in posSamples])
    return(scoreMap)

def scoreVariance(matData, matSamples):
    scoreMap = {}
    for feature in matData.keys():
        scoreMap[feature] = mData.mean_std([matData[feature][i] for i in matSamples], sample = True)[1]
    return(scoreMap)

def selectMutationNeighborhood(focusGene, mutSamples, dataFile, gPathway, trainSamples = None, tCutoff = 2.0, tIncrement = 0.5 , maxDistance = 2, method = "vsZero", penalty = 0.5):
    """performs simple univariate approach for selecting features with a threshold"""
    ## load data
    (matData, matFeatures, matSamples) = mData.rCRSData(dataFile, retFeatures = True)
    if trainSamples is not None and method != "variance":
        matSamples = list(set(matSamples) & set(trainSamples))
    mutatedSamples = list(set(mutSamples) & set(matSamples))
    nonSamples = list(set(matSamples) - set(mutSamples))
    permuteSamples = ["perm_%s" % (i) for i in mutatedSamples]
    
    ## score features by supervision
    nodeUnsigned = {}
    nodeScore = {}
    if method == "vsNon_svm":
        nodeUnsigned = scoreSVM(matData, mutatedSamples, nonSamples)
        nodeScore = rankScores(nodeUnsigned)
    if method == "vsPerm_svm":
        nodeUnsigned = scoreSVM(matData, mutatedSamples, permuteSamples)
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsNon":
        nodeUnsigned = scoreTT(matData, mutatedSamples, nonSamples)
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsPerm":
        nodeUnsigned = scoreTT(matData, mutatedSamples, permuteSamples)
        nodeScore = rankScores(nodeUnsigned)
    elif method == "variance":
        nodeUnsigned = scoreVariance(matData, matSamples)
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsMin":
        nonScores = rankScores(scoreTT(matData, mutatedSamples, nonSamples))
        varScores = rankScores(scoreVariance(matData, matSamples))
        for feature in list(set(nonScores.keys()) & set(varScores.keys())):
            minVal = 1.0
            valList = [nonScores[feature], varScores[feature]]
            for i in valList:
                try:
                    if abs(i) < abs(minVal):
                        minVal = i
                except TypeError:
                    pass
            nodeUnsigned[feature] = minVal
            nodeScore[feature] = minVal
    elif method == "vsMin_svm":
        nonScores = rankScores(scoreSVM(matData, mutatedSamples, nonSamples))
        varScores = rankScores(scoreVariance(matData, matSamples))
        for feature in list(set(nonScores.keys()) & set(varScores.keys())):
            minVal = 1.0
            valList = [nonScores[feature], varScores[feature]]
            for i in valList:
                try:
                    if abs(i) < abs(minVal):
                        minVal = i
                except TypeError:
                    pass
            nodeUnsigned[feature] = minVal
            nodeScore[feature] = minVal
    elif method == "vsMax":
        nonScores = rankScores(scoreTT(matData, mutatedSamples, nonSamples))
        varScores = rankScores(scoreVariance(matData, matSamples))
        for feature in list(set(nonScores.keys()) & set(varScores.keys())):
            maxVal = 0.0
            valList = [nonScores[feature], varScores[feature]]
            for i in valList:
                try:
                    if abs(i) > abs(maxVal):
                        maxVal = i
                except TypeError:
                    pass
            nodeUnsigned[feature] = maxVal
            nodeScore[feature] = maxVal
    elif method == "vsMax_svm":
        nonScores = rankScores(scoreSVM(matData, mutatedSamples, nonSamples))
        varScores = rankScores(scoreVariance(matData, matSamples))
        for feature in list(set(nonScores.keys()) & set(varScores.keys())):
            maxVal = 0.0
            valList = [nonScores[feature], varScores[feature]] 
            for i in valList:
                try:
                    if abs(i) > abs(maxVal):
                        maxVal = i
                except TypeError:
                    pass
                nodeUnsigned[feature] = maxVal
                nodeScore[feature] = maxVal
    
    ## initialize up and down pathways
    rinteractions = reverseInteractions(gPathway.interactions)
    upPathway = Pathway({focusGene : gPathway.nodes[focusGene]}, {})
    downPathway = Pathway({focusGene : gPathway.nodes[focusGene]}, {})
    coreNodes = set([focusGene])
    ## include complexes of focusGene
    if focusGene in gPathway.interactions:
        for target in gPathway.interactions[focusGene].keys():
            if gPathway.interactions[focusGene][target] == "component>":
                ## add complex to upstream pathway
                upPathway.nodes[target] = gPathway.nodes[target]
                if focusGene not in upPathway.interactions:
                    upPathway.interactions[focusGene] = {}
                upPathway.interactions[focusGene][target] = gPathway.interactions[focusGene][target]
                ## add complex to downstream pathway
                downPathway.nodes[target] = gPathway.nodes[target]
                if focusGene not in downPathway.interactions:
                    downPathway.interactions[focusGene] = {}
                downPathway.interactions[focusGene][target] = gPathway.interactions[focusGene][target]
                coreNodes.update([target])
    wPathway("core_pathway.tab", upPathway.nodes, upPathway.interactions)
    
    ## get frontier nodes up and down
    frontierUp = set()
    frontierDown = set()
    for node in upPathway.nodes.keys():
        frontierUp = frontierUp | set(getUpstream(node, maxDistance, gPathway.interactions))
        frontierDown = frontierDown | set(getDownstream(node, maxDistance, gPathway.interactions))
    frontierUp = frontierUp - frontierDown ##
    frontierUp = list(frontierUp - set([focusGene]))
    frontierDown = list(frontierDown - set([focusGene]))
    f = open("supervised.score", "w")
    for node in frontierUp:
        if node not in nodeScore:
            nodeScore[node] = -1
            nodeUnsigned[node] = "NA"
        f.write("%s\t%s\t%s\t+\n" % (node, nodeScore[node], nodeUnsigned[node]))
    for node in frontierDown:
        if node not in nodeScore:
            nodeScore[node] = -1
            nodeUnsigned[node] = "NA"
        f.write("%s\t%s\t%s\t-\n" % (node, nodeScore[node], nodeUnsigned[node]))
    f.close()
    frontierUp.sort(lambda x, y: cmp(nodeScore[y],nodeScore[x]))
    frontierDown.sort(lambda x, y: cmp(nodeScore[y],nodeScore[x]))
    
    ## add upPathway features by score
    count = 0
    hasUpstream = False
    tThresh = tCutoff
    for node in frontierUp:
        if ((count < 4) and (nodeScore[node] > 0)) or (nodeScore[node] > tThresh):
            addPaths = shortestPath(node, focusGene, gPathway.interactions)
            for currPath in addPaths:
                currPath.reverse()
                for addNode in currPath:
                    if addNode in upPathway.nodes:
                        continue
                    upPathway.nodes[addNode] = gPathway.nodes[addNode]
                    upPathway.interactions[addNode] = {}
                    if addNode in gPathway.interactions:
                        for target in gPathway.interactions[addNode].keys():
                            if target in upPathway.nodes:
                                upPathway.interactions[addNode][target] = gPathway.interactions[addNode][target] 
                    if addNode in rinteractions:
                        for source in rinteractions[addNode].keys():
                            if source in upPathway.nodes and source not in coreNodes:
                                upPathway.interactions[source][addNode] = gPathway.interactions[source][addNode] 
                    if addNode in matFeatures:
                        tThresh += tIncrement
                        count += 1
    if count >= 4:
        hasUpstream = True
    
    ## add downPathway features by score
    count = 0
    hasDownstream = False
    tThresh = tCutoff
    for node in frontierDown:
        if ((count < 4) and (nodeScore[node] > 0)) or (nodeScore[node] > tThresh):
            addPaths = shortestPath(focusGene, node, gPathway.interactions)
            for currPath in addPaths:
                for addNode in currPath:
                    if addNode in downPathway.nodes:
                        continue
                    downPathway.nodes[addNode] = gPathway.nodes[addNode]
                    if addNode in rinteractions:
                        for source in rinteractions[addNode].keys():
                            if source in downPathway.nodes:
                                if source not in downPathway.interactions:
                                    downPathway.interactions[source] = {}
                                downPathway.interactions[source][addNode] = gPathway.interactions[source][addNode]
                                if gPathway.interactions[source][addNode].startswith("-t"):
                                    tThresh += tIncrement
                                    count += 1
                    if addNode in gPathway.interactions:
                        for target in gPathway.interactions[addNode].keys():
                            if target in downPathway.nodes and target not in coreNodes:
                                if addNode not in downPathway.interactions:
                                    downPathway.interactions[addNode] = {}
                                downPathway.interactions[addNode][target] = gPathway.interactions[addNode][target] 
    if count >= 4:
        hasDownstream = True
    
    ## output pathway
    wPathway("upstream_pathway.tab", upPathway.nodes, upPathway.interactions)
    wPathway("downstream_pathway.tab", downPathway.nodes, downPathway.interactions)
    return(upPathway, downPathway, hasUpstream and hasDownstream)

def shiftCV(mutatedGene, mutatedSamples, dataSamples, trainSamples, uPathway, dPathway, 
            nNulls = 10, msepMethod = "tt", alpha = 0.05, mutsigFile = None):
    ## read paradigm output
    upScore = mData.rPARADIGM("%s_upstream.fa" % (mutatedGene), useRows = uPathway.nodes.keys())[1]
    downScore = mData.rPARADIGM("%s_downstream.fa" % (mutatedGene), useRows = dPathway.nodes.keys())[1]
    mData.wCRSData("paradigm_up.tab", upScore)
    mData.wCRSData("paradigm_down.tab", downScore)
    n_upScore = {}
    n_downScore = {}
    for null in range(1, nNulls+1):
        n_upScore[null] = mData.rPARADIGM("N%s_%s_upstream.fa" % (null, mutatedGene), useRows = [mutatedGene])[1]
        n_downScore[null] = mData.rPARADIGM("N%s_%s_downstream.fa" % (null, mutatedGene), useRows = [mutatedGene])[1]
    f = open("up.features", "w")
    for feature in (set(upScore[upScore.keys()[0]].keys()) - set([mutatedGene])):
        f.write("%s\n" % (feature))
    f.close()
    f = open("down.features", "w")
    for feature in (set(downScore[downScore.keys()[0]].keys()) - set([mutatedGene])):
        f.write("%s\n" % (feature))
    f.close()
    f = open("mut.features", "w")
    f.write("%s\n" % (mutatedGene))
    f.close()
    
    ## define sample groups
    trainGroup = trainSamples
    testGroup = list(set(dataSamples) - set(trainSamples))
    mutGroup = mutatedSamples
    nonGroup = list(set(dataSamples) - set(mutatedSamples))
    permGroup = ["perm_%s" % (sample) for sample in mutGroup]
    mutGroup_tr = list(set(mutGroup) & set(trainGroup))
    mutGroup_te = list(set(mutGroup) & set(testGroup))
    nonGroup_tr = list(set(nonGroup) & set(trainGroup))
    nonGroup_te = list(set(nonGroup) & set(testGroup))
    permGroup_tr = ["perm_%s" % (sample) for sample in mutGroup_tr]
    permGroup_te = ["perm_%s" % (sample) for sample in mutGroup_te]
    backGroup = []
    for sample in mutGroup:
        for null in range(1, nNulls+1):
            backGroup.append("null%s_%s" % (null, sample))
    f = open("include.samples", "w")
    for sample in (mutGroup_tr + nonGroup_tr):
        f.write("%s\n" % (sample))
    f.close()
    f = open("mut.circle", "w")
    f.write("id\t%s\n" % ("\t".join(mutGroup_tr + nonGroup_tr + permGroup_tr)))
    f.write("*")
    for sample in (mutGroup_tr + nonGroup_tr + permGroup_tr):
        if sample in mutGroup:
            f.write("\t1")
        elif sample in nonGroup:
            f.write("\t0")
        else:
            f.write("\t0.5")
    f.write("\n")
    f.close()
    
    ## score shifts      
    shiftScore = {}
    for sample in (mutGroup + nonGroup + permGroup):
        shiftScore[sample] = downScore[sample][mutatedGene] - upScore[sample][mutatedGene]
    for sample in (mutGroup + nonGroup + permGroup):
        for null in range(1, nNulls+1):
            shiftScore["null%s_%s" % (null, sample)] = n_downScore[null][sample][mutatedGene] - n_upScore[null][sample][mutatedGene]   
    f = open("shift.circle", "w")
    f.write("id\t%s\n" % ("\t".join(mutGroup_tr + nonGroup_tr + permGroup_tr)))
    f.write("*")
    for sample in (mutGroup_tr + nonGroup_tr + permGroup_tr):
        f.write("\t%s" % (shiftScore[sample]))
    f.write("\n")
    f.close()
    
    ## compute CV statistics
    (nonMean, nonStd) = mCalculate.mean_std([shiftScore[sample] for sample in nonGroup_tr])
    (permMean, permStd) = mCalculate.mean_std([shiftScore[sample] for sample in permGroup_tr])
    labelMap = {}
    nonScore_tr = {}
    nonScore_te = {}
    o = open("pshift.train.tab", "w")
    o.write("> %s\tP-Shifts:Table\n" % (mutatedGene))
    o.write("# sample\tclass\tP-Shift\n")
    f = open("mut.train.scores", "w")
    for sample in mutGroup_tr:
        labelMap[sample] = 1
        nonScore_tr[sample] = (shiftScore[sample]-nonMean)
        o.write("%s\t+\t%s\n" % (sample, nonScore_tr[sample]))
        f.write("%s\t%s\n" % (sample, nonScore_tr[sample]))
    f.close()
    f = open("non.train.scores", "w")
    for sample in nonGroup_tr:
        labelMap[sample] = 0
        nonScore_tr[sample] = (shiftScore[sample]-nonMean)
        o.write("%s\t-\t%s\n" % (sample, nonScore_tr[sample]))
        f.write("%s\t%s\n" % (sample, nonScore_tr[sample]))
    f.close()
    o.close()
    o = open("pshift.test.tab", "w")
    o.write("> %s\tP-Shifts:Table\n" % (mutatedGene))
    o.write("# sample\tclass\tP-Shift\n")
    f = open("mut.test.scores", "w")
    for sample in mutGroup_te:
        labelMap[sample] = 1
        nonScore_te[sample] = (shiftScore[sample]-nonMean)
        o.write("%s\t+\t%s\n" % (sample, nonScore_te[sample]))
        f.write("%s\t%s\n" % (sample, nonScore_te[sample]))
    f.close()
    f = open("non.test.scores", "w")
    for sample in nonGroup_te:
        labelMap[sample] = 0
        nonScore_te[sample] = (shiftScore[sample]-nonMean)
        o.write("%s\t-\t%s\n" % (sample, nonScore_te[sample]))
        f.write("%s\t%s\n" % (sample, nonScore_te[sample]))
    f.close()
    o.close()
    
    ## compute auc from stats
    sorted_tr = nonScore_tr.keys()
    sorted_tr.sort(lambda x, y: cmp(nonScore_tr[y],nonScore_tr[x]))
    auc_tr = computeAUC(sorted_tr, labelMap, nonScore_tr)[0]
    if len(testGroup) > 0:
        sorted_te = nonScore_te.keys()
        sorted_te.sort(lambda x, y: cmp(nonScore_te[y],nonScore_te[x]))
        auc_te = computeAUC(sorted_te, labelMap, nonScore_te)[0]
    else:
        auc_te = "NA"
    f = open("auc.stat", "w")
    f.write("%s\t%s\n" % (auc_tr, auc_te))
    f.close()
    
    ## compute mutant separation
    msepScore = {}
    msepScore["real"] = mutSeparation(nonGroup_tr, mutGroup_tr, shiftScore, method = msepMethod)
    for null in range(1, nNulls+1):
        msepScore["null%s" % (null)] = mutSeparation(["null%s_%s" % (null, sample) 
                                                     for sample in nonGroup_tr],
                                                     ["null%s_%s" % (null, sample)
                                                     for sample in mutGroup_tr],
                                                     shiftScore, method = msepMethod)
    f = open("real.scores", "w")
    f.write("%s\n" % (msepScore["real"]))
    f.close()
    f = open("null.scores", "w")
    for null in range(1, nNulls+1):
        f.write("%s\n" % (msepScore["null%s" % (null)]))
    f.close()
    
    ## compute background significance
    realScores = [msepScore["real"]]
    nullScores = [msepScore[null] for null in list(set(msepScore.keys()) - set(["real"]))]
    (nullMean, nullStd) = mCalculate.mean_std(nullScores)
    significanceScore = (realScores[0]-nullMean)/(nullStd+alpha)
    
    ## output sig.stats
    mutsigData = {}
    if mutsigFile is not None:
        mutsigData = mData.rCRSData(mutsigFile)["p"]
    if mutatedGene in mutsigData:
        mutVal = mutsigData[mutatedGene]
    else:
        mutVal = "NA"
    f = open("sig.tab", "w")
    f.write("# gene\tnon_mut\tmut\tmsep\tsignificance\tmutSig\n")
    f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (mutatedGene, len(nonGroup_tr), len(mutGroup_tr), 
                                      msepScore["real"], significanceScore, 
                                      str(mutVal).lstrip("<")))
    f.close()
    
def computeAUC(features, labelMap, scoreMap):
    points = []
    auc = 0.0
    x = 0.0
    index = 0
    while (index <= len(features)):
        (tp, fn, fp, tn) = getConfusion(index, features, labelMap)
        tpr = tp/float(tp+fn)
        fpr = fp/float(fp+tn)
        points.append( (fpr, tpr) )
        if fpr > x:
            dx = fpr - x
            x = fpr
            auc += dx*tpr
        index += 1
        if index < len(features):
            while (scoreMap[features[index-1]] == scoreMap[features[index]]):
                index += 1
                if index == len(features):
                    break
        elif index == len(features):
            pass
        else:
            break
    return(auc, points)

def getConfusion(index, features, labelMap):
    (tp, fn, fp, tn) = (0, 0, 0, 0)
    for feature in features[:index]:
        if labelMap[feature] == 1:
            tp += 1
        elif labelMap[feature] == 0:
            fp += 1
    for feature in features[index:]:
        if labelMap[feature] == 1:
            fn += 1
        elif labelMap[feature] == 0:
            tn += 1
    return(tp, fn, fp, tn)
