import math, os, sys, random, re
from copy import deepcopy
from array import array

def log(msg, file = None, die = False):
    """logger function"""
    if file is None:
        sys.stderr.write(msg)
    else:
        l = open(file, "a")
        l.write(msg)
        l.close()
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
        hugoList = rList(hugof)
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
    readPathway = Pathway({}, {})
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
                    readPathway.interactions[pline[1]] = {}
                if pline[0] not in readPathway.interactions[pline[1]]:
                    readPathway.interactions[pline[1]][pline[0]] = pline[2]
                else:
                    readPathway.interactions[pline[1]][pline[0]] += ";"+pline[2]
            else:
                if pline[0] not in readPathway.interactions:
                    readPathway.interactions[pline[0]] = {}
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
    componentMap = {}
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
                        outPathway.interactions[i] = {}
                    outPathway.interactions[i][j] = refInteractions[i][j]
    return(outPathway.nodes, outPathway.interactions)

def combinePathways(currentPathway, appendPathway, exclude = []):
    """combine Pathways"""
    combinedPathway = deepcopy(currentPathway)
    for source in appendPathway.interactions.keys():
        if source in exclude:
            continue
        if source not in combinedPathway.nodes:
            combinedPathway.nodes[source] = appendPathway.nodes[source]
        for target in appendPathway.interactions[source].keys():
            if target in exclude:
                continue
            if target not in combinedPathway.nodes:
                combinedPathway.nodes[target] = appendPathway.nodes[target]
            if source not in combinedPathway.interactions:
                combinedPathway.interactions[source] = {}
            if target not in combinedPathway.interactions[source]:
                combinedPathway.interactions[source][target] = appendPathway.interactions[source][target]
    return(combinedPathway)

def sortConnected(inPathway):
    rinteractions = reverseInteractions(inPathway.interactions)
    index = 1
    mapNets = {}
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
    netScore = {}
    for i in indexList:
        netScore[i] = len(mapNets[i])
    indexList.sort(lambda x, y: cmp(netScore[y], netScore[x]))
    for i in indexList:
        sortedNets.append(mapNets[i])
    return(sortedNets)

def flattenPathway(inPathway):
    """expands complexes into their respective gene components"""
    rinteractions = reverseInteractions(inPathway.interactions)
    allowedNodes = ["abstract", "miRNA", "protein", "rna"]
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
                    if inPathway.interactions[source][target] not in ["component>"]:
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
                                outPathway.interactions[source][element] = "-ca>"
                            else:
                                outPathway.interactions[source][element] = inPathway.interactions[source][target]
                elif inPathway.nodes[target] == "family":
                    if inPathway.interactions[source][target] not in ["component>", "member>"]:
                        #### some families have 0 members or complex members
                        for element in rinteractions[target].keys():
                            if rinteractions[target][element] == "member>":
                                if source not in outPathway.nodes:
                                    outPathway.nodes[source] = inPathway.nodes[source]
                                if element not in outPathway.nodes:
                                    outPathway.nodes[element] = inPathway.nodes[element]
                                if source not in outPathway.interactions:
                                    outPathway.interactions[source] = {}
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
            elif inPathway.nodes[source] == "family":
                if inPathway.nodes[target] in allowedNodes:
                    if inPathway.interactions[source][target] not in ["component>"]:
                        if source in rinteractions:
                            #### some families have 0 members or complex members
                            for element in rinteractions[source].keys():
                                if rinteractions[source][element] == "member>":
                                    if element not in outPathway.nodes:
                                        outPathway.nodes[element] = inPathway.nodes[element]
                                    if target not in outPathway.nodes:
                                        outPathway.nodes[target] = inPathway.nodes[target]
                                    if element not in outPathway.interactions:
                                        outPathway.interactions[element] = {}
                                    outPathway.interactions[element][target] = inPathway.interactions[source][target]
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

def signPath(path, interactions):
    sign = 1
    for i in range(len(path)-1):
        if interactions[path[i]][path[i+1]].endswith("|"):
            sign *= -1
    return(sign)

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

class DataMatrix:
    def __init__(self, inf, delim = "\t", null = "NA", useCols = None, useRows = None, 
                 type = "float"):
        ## initiate objects to be stored in class
        if type == "float":
            dataVals = array('f')
        else:
            dataVals = []
        colFeatures = []
        colIndex = {}
        rowFeatures = []
        rowIndex = {}
        
        ## read header
        f = open(inf, "r")
        line = f.readline().rstrip()
        if line.isspace():
            log("ERROR: encountered a blank on line 1\n", die = True)
        pline = line.split(delim)
        for ind, col in enumerate(pline[1:]):
            if useCols != None:
                if col not in useCols:
                    continue
            colIndex[col] = ind
            colFeatures.append(col)
        lineLength = len(pline)
        
        ## read data
        ind = 0
        for line in f:
            if line.isspace():
                continue
            line = line.rstrip()
            pline = line.split(delim)
            if len(pline) != lineLength:
                log("ERROR: length of line does not match the rest of %s\n" 
                    % (os.path.abspath(inf)), file = "err.log" , die = True)
            row = pline[0]
            if useRows != None:
                if row not in useRows:
                    continue
            rowIndex[row] = ind
            rowFeatures.append(row)
            for col in colFeatures:
                if (pline[colIndex[col]+1] == "") or (pline[colIndex[col]+1] == null):
                    if type == "float":
                        dataVals.append(float("nan"))
                    else:
                        dataVals.append(null)
                else: 
                    if type == "float":
                        dataVals.append(float(pline[colIndex[col]+1]))
                    else:
                        dataVals.append(pline[colIndex[col]+1])
            ind += 1
        f.close()
        self.values = deepcopy(dataVals)
        self.columns = deepcopy(colFeatures)
        self.cmap = deepcopy(colIndex)
        self.rows = deepcopy(rowFeatures)
        self.rmap = deepcopy(rowIndex)
    def get(self, col, row):
        rowIndex = self.rmap[row]
        colLength = len(self.columns)
        colIndex = self.cmap[col]
        return(self.values[rowIndex*colLength+colIndex])
        
def retColumns(inf, delim = "\t"):
    """returns the columns of a .tsv"""
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank header\n", die = True)
    line = line.rstrip("\r\n")
    return(re.split(delim, line)[1:])

def retRows(inf, delim = "\t", index = 0, header = True):
    """returns the rows of a .tsv"""
    rows = []
    f = open(inf, "r")
    if header:
        line = f.readline()
        if line.isspace():
            log("ERROR: encountered a blank header\n", die = True)
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        rows.append(re.split(delim, line)[index])
    return(rows)

def rCRSData(inf, appendData = {}, delim = "\t", null = "NA", useCols = None, useRows = None, retFeatures = False, debug = False):
    """reads .tsv into a [col][row] dictionary"""
    inData = {}  
    colFeatures = []
    rowFeatures = []
    ## copy appendData
    for col in appendData.keys():
        if useCols != None:
            if col not in useCols:
                continue
        inData[col] = {}
        for row in appendData[col].keys():
            if useRows != None:
                if row not in useRows:
                    continue
            inData[col][row] = [appendData[col][row]]
    ## read header
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    lineLength = len(pline)
    if debug:
        log("%s\nLENGTH: %s\n" % (line, lineLength))
    colIndex = {}
    for i, col in enumerate(pline[1:]):
        if useCols != None:
            if col not in useCols:
                continue
        colIndex[col] = i
        colFeatures.append(col)
        if col not in inData:
            inData[col] = {}
    ## read data
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        row = pline[0]
        if useRows != None:
            if row not in useRows:
                continue
        rowFeatures.append(row)
        if debug:
            log("%s\nLENGTH: %s\n" % (line, len(pline)))
        if len(pline) != lineLength:
            log("ERROR: length of line does not match the rest of %s\n" % (os.path.abspath(inf)), file = "err.log" , die = True)
        for col in colFeatures:
            if row not in inData[col]:
                inData[col][row] = []
            if pline[colIndex[col]+1] == "":
                inData[col][row].append(null)
            else:            
                inData[col][row].append(pline[colIndex[col]+1])
    f.close()
    ## average entries
    for col in inData.keys():
        for row in inData[col].keys():
            inData[col][row] = mean(inData[col][row], null = inData[col][row][0])
    if retFeatures:
        return(inData, colFeatures, rowFeatures)
    else:
        return(inData)

def wCRSData(outf, outData, delim = "\t", null = "NA", useCols = None, useRows = None):
    """write [col][row] dictionary to .tsv"""
    ## get colFeatures and rowFeatures
    if useCols == None:
        colFeatures = outData.keys()
    else:
        colFeatures = list(useCols)
    if useRows == None:
        rowFeatures = []
        for col in colFeatures:
            if col in outData:
                rowFeatures = outData[col].keys()
                break
    else:
        rowFeatures = list(useRows)
    ## write header
    if os.path.exists(outf):
        f = open(outf, "a")
    else:
        f = open(outf, "w")
        f.write("id%s\n" % (delim+delim.join(colFeatures)))
    ## write data
    for row in rowFeatures:
        f.write("%s" % (row))
        for col in colFeatures:
            try:
                f.write("%s" % (delim+str(outData[col][row])))
            except KeyError:
                f.write("%s" % (delim+null))
        f.write("\n")
    f.close()

def rwCRSData(outf, inf, delim = "\t", null = "NA", useCols = None, useRows = None, colMap = {}, rowMap = {}, rcMap = None):
    """reads and writes .tsv"""
    colFeatures = []
    ## read header
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
    lineLength = len(pline)
    colIndex = {}
    for i, col in enumerate(pline[1:]):
        colIndex[col] = i
        if useCols != None:
            if col not in useCols:
                continue
        colFeatures.append(col)
    ## write header
    if os.path.exists(outf):
        o = open(outf, "a")
    else:
        o = open(outf, "w")
        o.write("id%s\n" % (delim+delim.join(colFeatures)))
    ## read and write data
    rowCount = 0
    for line in f:
        if line.isspace():
            continue
        rowCount += 1
        line = line.rstrip("\r\n")
        pline = re.split(delim, line)
        if pline[0] in rowMap:
            mrow = rowMap[pline[0]]
        else:
            mrow = pline[0]
        if useRows != None:
            if mrow not in useRows:
                continue
        if len(pline) != lineLength:
            log("ERROR: length of line does not match the rest of the file\n", die = True)
        else:
            o.write("%s" % (mrow))
        if rcMap is not None:
            colMap = rcMap[mrow]
        for col in colFeatures:
            if col in colMap:
                mcol = colMap[col]
            else:
                mcol = col
            if pline[colIndex[mcol]+1] == "":
                o.write("%s" % (delim+null))
            else:            
                o.write("%s" % (delim+pline[colIndex[mcol]+1]))
        o.write("\n")
    f.close()
    o.close()

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

def rPARADIGM(inf, delim = "\t", useRows = None):
    """read PARADIGM format .fa output"""
    inLikelihood = {}
    inScore = {}
    f = open(inf, "r")
    for line in f:
        if line.isspace():
            continue
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            pline = re.split("[= ]", line)
            sample = pline[1]
            inLikelihood[sample] = float(pline[3])
            inScore[sample] = {}
        else:
            pline = re.split(delim, line)
            feature = pline[0]
            if useRows != None:
                if feature not in useRows:
                    continue
            inScore[sample][feature] = float(pline[1])
    f.close()
    return(inLikelihood, inScore)

def rMAF(inf, delim = "\t"):
    """read .maf format file"""
    mutData = {}
    mutClass = {}
    f = open(inf, "r")
    line = f.readline()
    if line.isspace():
        log("ERROR: encountered a blank on line 1\n", die = True)
    line = line.rstrip("\r\n")
    pline = re.split(delim, line)
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
        line = line.rstrip("\t\r\n")
        pline = re.split(delim, line)
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

def mean_std(inList, sample = True):
    """Calculates mean and std"""
    fList = floatList(inList)
    if len(fList) == 0:
        mean = "NA"
        std = "NA"
    else:
        mean = sum(fList)/float(len(fList))
        std = 0.0
        for i in fList:
            std += (i-mean)**2
        if len(fList) > 1:
            if sample:
                std = math.sqrt(std/(len(fList)-1))
            else:
                std = math.sqrt(std/len(fList))
        else:
            std = 0.0
    return(mean, std)

def getFeatureScores(dataMap, posSamples, negSamples, method = "variance",
                     outFile = None):
    """Computes feature scores for model selection"""
    ## define specific functions
    def scoreVariance(dataMatrix, useSamples):
        scoreMap = {}
        for feature in dataMatrix.columns:
            scoreMap[feature] = mean_std([dataMatrix.get(feature, i) for i in useSamples],
                                         sample = True)[1]
        return(scoreMap)
    def scoreTT(dataMatrix, posSamples, negSamples):
        scoreMap = {}
        for feature in dataMatrix.columns:
            scoreMap[feature] = ttest([dataMatrix.get(feature, i) for i in negSamples],
                                      [dataMatrix.get(feature, i) for i in posSamples])
        return(scoreMap)
    def scoreSVM(dataMatrix, posSamples, negSamples):
        import liblinear
        import liblinearutil
        svmLabels = []
        svmData = []
        featureMap = {}
        for sample in (posSamples + negSamples):
            if sample in posSamples:
                svmLabels.append(1)
            else:
                svmLabels.append(-1)
            svmData.append({})
            for i, feature in enumerate(dataMatrix.columns):
                if i+1 not in featureMap:
                    featureMap[i+1] = feature
                try:
                    svmData[-1][i+1] = dataMatrix.get(feature, sample)
                except ValueError:
                    svmData[-1][i+1] = 0.0
        prob = liblinear.problem(svmLabels, svmData)
        param = liblinear.parameter('-s 3 -c 5 -q')
        liblinearutil.save_model('model_file', liblinearutil.train(prob, param))
        # m = liblinearutil.train(prob, param)
        # testLabels = [] ## like svmLabels
        # testData = [] ## like svmData
        # for i in range(len(testLabels)):
        #     x0, max_idx = gen_feature_nodearray(testData[i])
        #     label = liblinear.predict(m, x0)
        #     if label == testLabels[i]:
        #         print "correct!"
        weights = rList("model_file")[6:]
        scoreMap = {}
        for feature in featureMap.keys():
            scoreMap[featureMap[feature]] = float(weights[feature-1])
        return (scoreMap)
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
        for attachment in scoreMap:
            rankMap[attachment] = {}
        for attachment in scoreMap:
            rankFeatures = []
            for feature in scoreMap[attachment]:
                try:
                    fval = float(scoreMap[attachment][feature])
                    if fval != fval:
                        raise ValueError
                    rankFeatures.append(feature)
                except ValueError:
                    rankMap[attachment][feature] = 0
            rankFeatures.sort(lambda x, y: cmp(abs(scoreMap[attachment][x]),abs(scoreMap[attachment][y])))
            for i, feature in enumerate(rankFeatures):
                rankMap[attachment][feature] = float(i+1)/len(rankFeatures)
        return(rankMap)
    
    ## score features by method
    varianceScore = {"mRNA" : {}, "active" : {}}
    ttestScore = {"mRNA" : {}, "active" : {}}
    svmScore = {"mRNA" : {}, "active" : {}}
    for attachment in dataMap:
        matData = DataMatrix(dataMap[attachment])
        if method in ["variance", "vsMin", "vsMin_svm", "vsMax", "vsMax_svm"]:
            varianceScore[attachment] = scoreVariance(matData, posSamples+negSamples)
        if method in ["vsNon", "vsMin", "vsMax"]:
            ttestScore[attachment] = scoreTT(matData, posSamples, negSamples)
        if method in ["vsNon_svm", "vsMin_svm", "vsMax_svm"]:
            svmScore[attachment] = scoreSVM(matData, posSamples, negSamples)
    
    ## compute score ranks
    nodeUnsigned = {}
    nodeScore = {}
    if method == "variance":
        nodeUnsigned = varianceScore
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsNon":
        nodeUnsigned = ttestScore
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsNon_svm":
        nodeUnsigned = svmScore
        nodeScore = rankScores(nodeUnsigned)
    elif method == "vsMin":
        nodeUnsigned = ttestScore
        nonScores = rankScores(ttestScore)
        varScores = rankScores(varianceScore)
        for attachment in nodeUnsigned:
            if attachment not in nodeScore:
                nodeScore[attachment] = {}
            for feature in nodeUnsigned[attachment]:
                minVal = 1.0
                valList = [nonScores[attachment][feature], varScores[attachment][feature]]
                for i in valList:
                    try:
                        if abs(i) < abs(minVal):
                            minVal = i
                    except TypeError:
                        pass
                nodeScore[attachment][feature] = minVal
    elif method == "vsMin_svm":
        nodeUnsigned = svmScore
        nonScores = rankScores(svmScore)
        varScores = rankScores(varianceScore)
        for attachment in nodeUnsigned:
            if attachment not in nodeScore:
                nodeScore[attachment] = {}
            for feature in nodeUnsigned[attachment]:
                minVal = 1.0
                valList = [nonScores[attachment][feature], varScores[attachment][feature]]
                for i in valList:
                    try:
                        if abs(i) < abs(minVal):
                            minVal = i
                    except TypeError:
                        pass
                nodeScore[attachment][feature] = minVal
    elif method == "vsMax":
        nodeUnsigned = ttestScore
        nonScores = rankScores(ttestScore)
        varScores = rankScores(varianceScore)
        for attachment in nodeUnsigned:
            if attachment not in nodeScore:
                nodeScore[attachment] = {}
            for feature in nodeUnsigned[attachment]:
                maxVal = 0.0
                valList = [nonScores[attachment][feature], varScores[attachment][feature]]
                for i in valList:
                    try:
                        if abs(i) > abs(maxVal):
                            maxVal = i
                    except TypeError:
                        pass
                nodeScore[attachment][feature] = maxVal
    elif method == "vsMax_svm":
        nodeUnsigned = svmScore
        nonScores = rankScores(svmScore)
        varScores = rankScores(varianceScore)
        for attachment in nodeUnsigned:
            if attachment not in nodeScore:
                nodeScore[attachment] = {}
            for feature in nodeUnsigned[attachment]:
                maxVal = 0.0
                valList = [nonScores[attachment][feature], varScores[attachment][feature]]
                for i in valList:
                    try:
                        if abs(i) > abs(maxVal):
                            maxVal = i
                    except TypeError:
                        pass
                nodeScore[attachment][feature] = maxVal
    
    ## output scores
    if outFile:
        o = open(outFile, "w")
        o.write("id\tattachment\tvalue\trank\n")
        for attachment in nodeScore:
            for feature in nodeScore[attachment]:
                o.write("%s\t%s\t%s\t%s\n" % (feature, attachment, nodeUnsigned[attachment][feature], nodeScore[attachment][feature]))
        o.close()
    return(nodeUnsigned, nodeScore)

def identifyNeighborhoods(focusGene, globalPathway, maxDistance = 2):
    """this"""
    gNodes = globalPathway.nodes
    fgInteractions = globalPathway.interactions
    rgInteractions = reverseInteractions(fgInteractions)
    
    ## define specific functions
    def globComplexes(focusGene, gNodes, fgInteractions, rgInteractions):
        if focusGene not in fgInteractions:
            return([], {})
        complexes = []
        for target in fgInteractions[focusGene]:
            if fgInteractions[focusGene][target] == "component>":
                complexes.append(target)
        seenComplexes = set()
        groupedComplexes = []
        isDownstream = {}
        while len(complexes) > 0:
            currentComplex = complexes.pop(0)
            isDownstream[currentComplex] = True
            if currentComplex in seenComplexes:
                continue
            seenComplexes.update([currentComplex])
            currentGroup = [currentComplex]
            currentUnwalked = [currentComplex]
            while len(currentUnwalked) > 0:
                currentNode = currentUnwalked.pop(0)
                if currentNode in fgInteractions:
                    for target in fgInteractions[currentNode]:
                        if fgInteractions[currentNode][target] == "component>" and gNodes[target] == "complex" and target not in seenComplexes:
                            currentGroup.append(target)
                            currentUnwalked.append(target)
                            seenComplexes.update([target])
                if currentNode in rgInteractions:
                    for source in rgInteractions[currentNode]:
                        if rgInteractions[currentNode][source] == "component>" and gNodes[source] == "complex" and source not in seenComplexes:
                            currentGroup.append(source)
                            currentUnwalked.append(source)
                            seenComplexes.update([source])
            groupedComplexes.append(currentGroup)
        assignCount = 1
        while assignCount > 0:
            assignCount = 0
            assignFeatures = isDownstream.keys()
            for complex in assignFeatures:
                if complex in fgInteractions:
                    for target in fgInteractions[complex]:
                        if target in seenComplexes and target not in isDownstream:
                            isDownstream[target] = True
                            assignCount += 1
        for complex in seenComplexes:
            if complex not in isDownstream:
                isDownstream[complex] = False
        return(groupedComplexes, isDownstream)
    def searchUpstream(focusGene, focusFeatures, gNodes, rgInteractions, maxDistance = 2):
        upstreamPathway = Pathway({focusGene : gNodes[focusGene]}, {})
        upstreamFeatures = {"active" : []}
        seenNodes = set(focusFeatures)
        searchNodes = deepcopy(focusFeatures)
        frontierNodes = []
        for i in range(1, maxDistance+1):
            while(len(searchNodes) > 0):
                currentNode = searchNodes.pop(0)
                if currentNode in rgInteractions:
                    for source in rgInteractions[currentNode]:
                        if source == focusGene:
                            continue
                        elif i == 1:
                            if gNodes[source] == "abstract":
                                continue
                            elif rgInteractions[currentNode][source] == "component>":
                                upstreamPathway.nodes[source] = gNodes[source]
                                if source not in upstreamPathway.interactions:
                                    upstreamPathway.interactions[source] = {}
                                upstreamPathway.interactions[source][focusGene] = "-a>"
                                upstreamFeatures["active"].append(source)
                                seenNodes.update([source])
                                frontierNodes.append(source)
                            elif rgInteractions[currentNode][source].startswith("-a"):
                                upstreamPathway.nodes[source] = gNodes[source]
                                if source not in upstreamPathway.interactions:
                                    upstreamPathway.interactions[source] = {}
                                upstreamPathway.interactions[source][focusGene] = rgInteractions[currentNode][source]
                                upstreamFeatures["active"].append(source)
                                seenNodes.update([source])
                                frontierNodes.append(source)
                            elif rgInteractions[currentNode][source].startswith("-t"):
                                upstreamPathway.nodes[source] = gNodes[source]
                                if source not in upstreamPathway.interactions:
                                    upstreamPathway.interactions[source] = {}
                                upstreamPathway.interactions[source][focusGene] = rgInteractions[currentNode][source]
                                seenNodes.update([source])
                        elif source in upstreamPathway.nodes:
                            upstreamPathway.interactions[source][currentNode] = rgInteractions[currentNode][source]
                        else:
                            if gNodes[source] == "abstract":
                                continue
                            else:
                                upstreamPathway.nodes[source] = gNodes[source]
                                if source not in upstreamPathway.interactions:
                                    upstreamPathway.interactions[source] = {}
                                upstreamPathway.interactions[source][currentNode] = rgInteractions[currentNode][source]
                                upstreamFeatures["active"].append(source)
                                seenNodes.update([source])
                                frontierNodes.append(source)
            searchNodes = deepcopy(frontierNodes)
        return(upstreamFeatures, upstreamPathway)
    def searchDownstream(focusGene, focusFeatures, gNodes, fgInteractions, maxDistance = 2):
        downstreamPathway = Pathway({focusGene : gNodes[focusGene]}, {})
        downstreamFeatures = {"mRNA" : [], "active" : []}
        seenNodes = set(focusFeatures)
        searchNodes = deepcopy(focusFeatures)
        frontierNodes = []
        for i in range(1, maxDistance+1):
            while(len(searchNodes) > 0):
                currentNode = searchNodes.pop(0)
                if currentNode in fgInteractions:
                    for target in fgInteractions[currentNode]:
                        if target == focusGene:
                            continue
                        elif i == 1:
                            if gNodes[target] == "abstract":
                                continue
                            elif fgInteractions[currentNode][target] == "component>":
                                continue
                            elif fgInteractions[currentNode][target].startswith("-a"):
                                downstreamPathway.nodes[target] = gNodes[target]
                                if focusGene not in downstreamPathway.interactions:
                                    downstreamPathway.interactions[focusGene] = {}
                                downstreamPathway.interactions[focusGene][target] = fgInteractions[currentNode][target]
                                downstreamFeatures["active"].append(target)
                                seenNodes.update([target])
                                frontierNodes.append(target)
                            elif fgInteractions[currentNode][target].startswith("-t"):
                                downstreamPathway.nodes[target] = gNodes[target]
                                if focusGene not in downstreamPathway.interactions:
                                    downstreamPathway.interactions[focusGene] = {}
                                downstreamPathway.interactions[focusGene][target] = fgInteractions[currentNode][target]
                                downstreamFeatures["mRNA"].append(target)
                                seenNodes.update([target])
                        elif target in downstreamPathway.nodes:
                            if currentNode not in downstreamPathway.interactions:
                                downstreamPathway.interactions[currentNode] = {}
                            downstreamPathway.interactions[currentNode][target] = fgInteractions[currentNode][target]
                        else:
                            if gNodes[target] == "abstract":
                                continue
                            else:
                                downstreamPathway.nodes[target] = gNodes[target]
                                if currentNode not in downstreamPathway.interactions:
                                    downstreamPathway.interactions[currentNode] = {}
                                downstreamPathway.interactions[currentNode][target] = fgInteractions[currentNode][target]
                                if fgInteractions[currentNode][target].startswith("-t"):
                                    downstreamFeatures["mRNA"].append(target)
                                else:
                                    frontierNodes.append(target)
                                seenNodes.update([target])
            searchNodes = deepcopy(frontierNodes)
        return(downstreamFeatures, downstreamPathway)
    
    ## check if focusGene belongs to any complexes
    (complexGroups, complexDownstream) = globComplexes(focusGene, gNodes, fgInteractions, rgInteractions)
    
    ## identify upstream and downstream of focusGene and complexes
    focusPathways = []
    (upstreamFeatures, upstreamPathway) = searchUpstream(focusGene, [focusGene], gNodes, rgInteractions, maxDistance = maxDistance)
    (downstreamFeatures, downstreamPathway) = searchDownstream(focusGene, [focusGene], gNodes, fgInteractions, maxDistance = maxDistance)
    focusPathways.append([[upstreamFeatures, upstreamPathway], [downstreamFeatures, downstreamPathway]])
    for currentGroup in complexGroups:
        currentDownstream = []
        for complex in currentGroup:
            if complexDownstream[complex]:
                currentDownstream.append(complex)
        (upstreamFeatures, upstreamPathway) = searchUpstream(focusGene, currentGroup, gNodes, rgInteractions, maxDistance = maxDistance)
        (downstreamFeatures, downstreamPathway) = searchDownstream(focusGene, currentDownstream, gNodes, fgInteractions, maxDistance = maxDistance)
        focusPathways.append([[upstreamFeatures, upstreamPathway], [downstreamFeatures, downstreamPathway]])
    return(focusPathways)

def selectMutNeighborhood(focusGene, nodeScore, nodeRank, upstreamFeatures, downstreamFeatures, upstreamBase, downstreamBase, threshold = 0.68, penalty = 0.01):
    """performs simple univariate approach for selecting features with a threshold"""
    upstreamPathway = Pathway({focusGene : upstreamBase.nodes[focusGene]}, {})
    addedFeatures = []
    log("## upstream\n", file = "selection.log")
    rankMap = {}
    rankFeatures = []
    for feature in upstreamFeatures["active"]:
        if feature in nodeRank["active"]:
            rankMap[feature] = nodeRank["active"][feature]
            if feature not in rankFeatures:
                rankFeatures.append(feature)
        elif feature in nodeRank["mRNA"]:
            rankMap[feature] = nodeRank["mRNA"][feature]
            if feature not in rankFeatures:
                rankFeatures.append(feature)
    rankFeatures.sort(lambda x, y: cmp(rankMap[y], rankMap[x]))
    while (len(rankFeatures) > 0) and (len(addedFeatures) < 4 or rankMap[rankFeatures[0]] >= threshold+penalty*len(addedFeatures)):
        currentFeature = rankFeatures.pop(0)
        addedFeatures.append(currentFeature)
        addPaths = shortestPath(currentFeature, focusGene, upstreamBase.interactions)
        for path in addPaths:
            for i in range(len(path)-1):
                upstreamPathway.nodes[path[i]] = upstreamBase.nodes[path[i]]
                if path[i] not in upstreamPathway.interactions:
                    upstreamPathway.interactions[path[i]] = {}
                upstreamPathway.interactions[path[i]][path[i+1]] = upstreamBase.interactions[path[i]][path[i+1]]
    for feature in addedFeatures + rankFeatures:
        featureVals = []
        if feature in nodeRank["mRNA"]:
            featureVals.append("mRNA:%s(%s)" % (nodeScore["mRNA"][feature], nodeRank["mRNA"][feature]))
        if feature in nodeRank["active"]:
            featureVals.append("active:%s(%s)" % (nodeScore["active"][feature], nodeRank["active"][feature]))
        if feature in addedFeatures:
            log("%s\t%s\t%s\n" % (feature, ",".join(featureVals), "*"), file = "selection.log")
        elif feature in upstreamPathway.nodes:
            log("%s\t%s\t%s\n" % (feature, ",".join(featureVals), "-"), file = "selection.log")
        else:
            log("%s\t%s\n" % (feature, ",".join(featureVals)), file = "selection.log")
    downstreamPathway = Pathway({focusGene : downstreamBase.nodes[focusGene]}, {})
    addedFeatures = []
    log("## downstream\n", file = "selection.log")
    rankMap = {}
    rankFeatures = []
    for feature in downstreamFeatures["active"]:
        if feature in nodeRank["active"]:
            rankMap[feature] = nodeRank["active"][feature]
            if feature not in rankFeatures:
                rankFeatures.append(feature)
    for feature in downstreamFeatures["mRNA"]:
        if feature in nodeRank["mRNA"]:
            rankMap[feature] = nodeRank["mRNA"][feature]
            if feature not in rankFeatures:
                rankFeatures.append(feature)
    rankFeatures.sort(lambda x, y: cmp(rankMap[y],rankMap[x]))
    while (len(rankFeatures) > 0) and (len(addedFeatures) < 4 or rankMap[rankFeatures[0]] >= threshold+penalty*len(addedFeatures)):
        currentFeature = rankFeatures.pop(0)
        addedFeatures.append(currentFeature)
        addPaths = shortestPath(focusGene, currentFeature, downstreamBase.interactions)
        for path in addPaths:
            for i in range(len(path)-1):
                downstreamPathway.nodes[path[i+1]] = downstreamBase.nodes[path[i+1]]
                if path[i] not in downstreamPathway.interactions:
                    downstreamPathway.interactions[path[i]] = {}
                downstreamPathway.interactions[path[i]][path[i+1]] = downstreamBase.interactions[path[i]][path[i+1]]
    for feature in addedFeatures + rankFeatures:
        featureVals = []
        if feature in nodeRank["mRNA"]:
            featureVals.append("mRNA:%s(%s)" % (nodeScore["mRNA"][feature], nodeRank["mRNA"][feature]))
        if feature in nodeRank["active"]:
            featureVals.append("active:%s(%s)" % (nodeScore["active"][feature], nodeRank["active"][feature]))
        if feature in addedFeatures:
            log("%s\t%s\t%s\n" % (feature, ",".join(featureVals), "*"), file = "selection.log")
        elif feature in downstreamPathway.nodes:
            log("%s\t%s\t%s\n" % (feature, ",".join(featureVals), "-"), file = "selection.log")
        else:
            log("%s\t%s\n" % (feature, ",".join(featureVals)), file = "selection.log")
    return(upstreamPathway, downstreamPathway)

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

def computeAUC(features, labelMap, scoreMap):
    points = []
    auc = 0.0
    x = 0.0
    index = 0
    while (index <= len(features)):
        (tp, fn, fp, tn) = getConfusion(index, features, labelMap)
        try:
            tpr = tp/float(tp+fn)
            fpr = fp/float(fp+tn)
        except ZeroDivisionError:
            return("NA", [])
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

def ttest(values0, values1, alpha = 0.05, retDegrees = False):
    (mean0, std0) = mean_std(values0)
    (mean1, std1) = mean_std(values1)
    try:
        tval = (mean1-mean0)/(math.sqrt((std1**2)/len(values1)+(std0**2)/len(values0))+alpha)
    except ZeroDivisionError:
        tval = "NA"
    except TypeError:
        tval = "NA"
    if retDegrees:
        df = (((std0**2)/len(values0)+(std1**2)/len(values1))**2)/((((std0**2)/len(values0))**2)/(len(values0)-1)+(((std1**2)/len(values1))**2)/(len(values1)-1))
        return(tval, df)
    else:
        return(tval)

def kstest(values0, values1):
    (mean0, std0) = mean_std(values0)
    (mean1, std1) = mean_std(values1)
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

