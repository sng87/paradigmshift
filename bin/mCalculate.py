## Calculate module
## Written By: Sam Ng
## Last Updated: 5/17/11
import math, sys
import mData

def log(msg, die = False):
    sys.stderr.write(msg)
    if die:
        sys.exit(1)

def revDict(inDict, list = True):
    outDict = dict()
    for i in inDict.keys():
        if inDict[i] in outDict:
            if list:
                outDict[inDict[i]].append(str(i))
            else:
                outDict[inDict[i]] = outDict[inDict[i]]+";"+str(i)
        else:
            if list:
                outDict[inDict[i]] = [str(i)]
            else:
                outDict[inDict[i]] = str(i)
    return(outDict)

def removeFeature(feature, inList):
    outList = []
    for i in inList:
        if feature == i:
            continue
        outList.append(i)
    return(outList)

def getValues(features, inDict):
    outList = []
    for i in features:
        if i not in inDict:
            continue
        outList.append(inDict[i])
    return(outList)

def proportions(inList):
    pMap = dict()
    for i in set(inList):
        pMap[i] = 0
    for i in inList:
        pMap[i] += 1
    for i in set(inList):
        pMap[i] = pMap[i]/float(len(inList))
    return(pMap)

def mean(inList, null = "NA"):
    """Calculates mean"""
    
    cList = mData.floatList(inList)
    if len(cList) == 0:
        mean = null
    else:
        mean = sum(cList)/len(cList)
    return (mean)

def mean_std(inList, sample = True):
    """Calculates mean and std"""
    
    cList = mData.floatList(inList)
    if len(cList) == 0:
        mean = "NA"
        std = "NA"
    else:
        mean = sum(cList)/float(len(cList))
        std = 0.0
        for i in cList:
            std += (i-mean)**2
        if len(cList) > 1:
            if sample:
                std = math.sqrt(std/(len(cList)-1))
            else:
                std = math.sqrt(std/len(cList))
        else:
            std = 0.0
    return(mean, std)

def median(inList):
    """Calculates median"""
    
    cList = mData.floatList(inList)
    cList.sort()
    if len(cList) == 0:
        median = "NA"
    else:
        if len(cList)%2 == 1:
            median = cList[len(cList)/2]
        else:
            median = (cList[len(cList)/2]+cList[(len(cList)/2)-1])/2.0
    return(median)

def quartiles(inList):
    """Returns the 25/50/75 quartiles"""
    
    cList = mData.floatList(inList)
    cList.sort()
    if len(cList) < 2:
        boundaries = ["NA", "NA", "NA"]
    else:
        boundaries = [median(cList[:len(cList)/2]), median(cList), median(cList[len(cList)/2:])]
    return(boundaries)

def pcorrelation(list1, list2):
    """Calculates pearson correlation"""
        
    if len(list1) != len(list2):
        log("ERROR: sizes of list are not equal\n", die = True)
    mean1 = mean(list1)
    mean2 = mean(list2)
    cov = 0.0
    stdev1 = 0.0
    stdev2 = 0.0
    for i in range(len(list1)):
        try:
            fval1 = float(list1[i])
            fval2 = float(list2[i])
            cov += (fval1-mean1)*(fval2-mean2)
            stdev1 += (fval1-mean1)**2
            stdev2 += (fval2-mean2)**2
        except ValueError:
            continue
    stdev1 = math.sqrt(stdev1)
    stdev2 = math.sqrt(stdev2)
    if stdev1 == 0 or stdev2 == 0:
        value = "NA"
    else:
        value = cov/(stdev1*stdev2)
    return(value)

def correlationPairwise(outf, inf, group1, group2, method = pcorrelation):
    """Computes all pairwise correlations between group1 and group2"""
    
    inData = mData.rCRSData(inf)
    inRows = inData[inData.keys()[0]].keys()
    f = open(outf, "w")
    for i in group1:
        if i not in inData:
            continue
        list1 = []
        for k in inRows:
            list1.append(inData[i][k])
        for j in group2:
            if j not in inData:
                continue
            list2 = []
            for k in inRows:
                list2.append(inData[j][k])
            value = method(list1, list2)
            f.write("%s\t%s\t%s\n" % (i, j, value))
    f.close()

def correlationMatrix(outf, inf, method = pcorrelation):
    """Takes a tab file and cross-correlates the columns"""
    
    outData = dict()
    inData = mData.rCRSData(inf)
    inCols = inData.keys()
    inRows = inData[inCols[0]].keys()
    for i in inCols:
        outData[i] = dict()
        outData[i][i] = 1.0
    for i in range(len(inCols)-1):
        list1 = []
        for k in inRows:
            list1.append(inData[inCols[i]][k])
        for j in range(i+1, len(inCols)):
            list2 = []
            for k in inRows:
                list2.append(inData[inCols[j]][k])
            value = method(list1, list2)
            outData[inCols[i]][inCols[j]] = value
            outData[inCols[j]][inCols[i]] = value
    mData.wCRSData(outf, outData)       

def ttest(group0, group1, inData, alpha = 0.01):
    """Computes t-statistics between two groups"""
    
    scoreMap = dict()
    for i in inData[inData.keys()[0]].keys():
        values0 = []
        values1 = []
        for j in inData.keys():
            if inData[j][i] == "NA":
                continue
            if j in group0:
                values0.append(inData[j][i])
            elif j in group1:
                values1.append(inData[j][i])
        if (len(values0) < 2) | (len(values1) < 2):
            scoreMap[i] = "NA"
        else:
            (mean0, std0) = mean_std(values0)
            (mean1, std1) = mean_std(values1)
            if (std0 == 0) & (std1 == 0):
                scoreMap[i] = "NA"
            else:
                scoreMap[i] = (mean0-mean1)/(math.sqrt((std0**2)/len(values0)+(std1**2)/len(values1))+alpha)
    return(scoreMap)

def zIPL(inData, method = "deviation"):
    """Computes z-scores for deviations from 0"""
    
    scoreMap = dict()
    for i in inData[inData.keys()[0]].keys():
        values = []
        for j in inData.keys():
            if inData[j][i] == "NA":
                continue
            else:
                values.append(inData[j][i])
        if len(values) < 2:
            scoreMap[i] = 0
        else:
            (mean, std) = mean_std(values)
            if std == 0:
                scoreMap[i] = 0
            else:
                scoreMap[i] = (mean)/(std)
    return(scoreMap)

def sign(value):
    """Returns the sign"""
    
    if value == 0:
        return(0)
    elif value > 0:
        return(1)
    elif value < 0:
        return(-1)

def computeROC(positivef, negativef):
    """Computes AUC and ROC curve"""
    
    allData = mData.r2Col(positivef)
    posSamples = set(allData.keys())
    allData = mData.r2Col(negativef, appendData = allData)
    negSamples = set(allData.keys()) - posSamples
    allSamples = allData.keys()
    allSamples.sort(lambda x, y: cmp(allData[y],allData[x]))
    x = [0.0]
    y = [0.0]
    TP = 0
    FP = 0
    TN = len(negSamples)
    FN = len(posSamples)
    for i in allSamples:
        if i in posSamples:
            TP += 1
            FN -= 1
        else:
            TN -= 1
            FP += 1
        x.append(float(FP)/len(negSamples))
        y.append(float(TP)/len(posSamples))
    AUC = 0.0
    for i in range(len(x)-1):
        AUC += (x[i+1]-x[i])*((y[i+1]+y[i])/2)
    return(AUC, x, y)
