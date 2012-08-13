#!/usr/bin/env python
"""layout-cytoscapeweb.py: Generate cytoscapeweb sessions from .sif and .NA files

Usage:
  layout-cytoscapeweb.py [options] feature destDir

Options:
  -n str    note text
  -s        output score table
  -q        run quietly
"""
## Written By: Sam Ng
## Last Updated: 5/18/11
import os, os.path, sys, getopt, re
import mData, mPathway, mCalculate

verbose = True
netExtension = ".sif"
noteText = "Jorma"

htmlLink = """<a href="%s">%s</a>"""

htmlHead = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html>
    
    <head>
        <title>CytoscapeWeb - UCSC PARADIGM Subnets</title>
        
        <script type="text/javascript" src="../js/min/json2.min.js"></script>
        <script type="text/javascript" src="../js/min/AC_OETags.min.js"></script>
        <script type="text/javascript" src="../js/min/cytoscapeweb.min.js"></script>
        
        <script type="text/javascript">
            window.onload = function() {
                // id of Cytoscape Web container div
                var div_id = "cytoscapeweb";
                
                // NOTE: - the attributes on nodes and edges
                //       - it also has directed edges, which will automatically display edge arrows
                var xml = '\\
                """
htmlTail = """                ';
                
                // visual style we will use
                var visual_style = {
                    global: {
                        backgroundColor: "#FFFFFF"
                    },
                    nodes: {
                        shape: {
                            discreteMapper: {
                                attrName: "type",
                                entries: [
                                    { attrValue: "protein", value: "CIRCLE" },
                                    { attrValue: "complex", value: "HEXAGON" },
                                    { attrValue: "abstract", value: "ROUNDRECT" },
                                    { attrValue: "drug", value: "TRIANGLE" },
                                    { attrValue: "plot", value: "CIRCLE" }
                                ]
                            }
                        },
                        borderWidth: 2,
                        borderColor: "#0000000",
                        size: {
                            passthroughMapper: { attrName: "size" }
                        },
                        color: {
                            passthroughMapper: { attrName: "color" }
                        },
                        image: {
                            passthroughMapper: { attrName: "image" }
                        },
                        labelHorizontalAnchor: "center"
                    },
                    edges: {
                        width: 2,
                        color: "#000000",
                        targetArrowShape: {
                            discreteMapper: {
                                attrName: "interaction",
                                    entries: [
                                    { attrValue: "-t|", value: "T" },
                                    { attrValue: "-t>", value: "ARROW" },
                                    { attrValue: "-a|", value: "T" },
                                    { attrValue: "-a>", value: "ARROW" },
                                    { attrValue: "-ap|", value: "T" },
                                    { attrValue: "-ap>", value: "ARROW" },
                                    { attrValue: "component>", value: "NONE" },
                                    { attrValue: "-disconnected-", value: "NONE" }
                                ]
                            }
                        },
                        style: {
                            discreteMapper: {
                                attrName: "interaction",
                                entries: [
                                    { attrValue: "-t|", value: "SOLID" },
                                    { attrValue: "-t>", value: "SOLID" },
                                    { attrValue: "-a|", value: "LONG_DASH" },
                                    { attrValue: "-a>", value: "LONG_DASH" },
                                    { attrValue: "-ap|", value: "LONG_DASH" },
                                    { attrValue: "-ap>", value: "LONG_DASH" },
                                    { attrValue: "component>", value: "LONG_DASH" },
                                    { attrValue: "-disconnected-", value: "LONG_DASH" }
                                ]
                            }
                        }
                    }
                };
                
                // initialization options
                var options = {
                    swfPath: "../swf/CytoscapeWeb",
                    flashInstallerPath: "../swf/playerProductInstall"
                };
                
                var vis = new org.cytoscapeweb.Visualization(div_id, options);
                
                vis.ready(function() {
                    // add a listener for when nodes and edges are clicked
                    vis.addListener("click", "nodes", function(event) {
                        handle_click(event);
                    })
                    .addListener("click", "edges", function(event) {
                        handle_click(event);
                    });
                    
                    function handle_click(event) {
                         var target = event.target;
                         
                         clear();
                         print("event.group = " + event.group);
                         for (var i in target.data) {
                            var variable_name = i;
                            var variable_value = target.data[i];
                            print( "event.target.data." + variable_name + " = " + variable_value );
                         }
                    }
                    
                    function clear() {
                        document.getElementById("note").innerHTML = "";
                    }
                    
                    function print(msg) {
                        document.getElementById("note").innerHTML += "<p>" + msg + "</p>";
                    }
                });

                var draw_options = {
                    // your data goes here
                    network: xml,
                    // hide edge labels
                    edgeLabelsVisible: false,
                    // let's try another layout
                    layout: "ForceDirected",
                    // set the style at initialisation
                    visualStyle: visual_style,
                    // show pan zoom
                    panZoomControlVisible: true 
                };
                
                vis.draw(draw_options);
            };
        </script>
        
        <style type= >
            * { margin: 0; padding: 0; font-family: Helvetica, Arial, Verdana, sans-serif; }
            html, body { height: 100%; width: 100%; padding: 0; margin: 0; background-color: #f0f0f0; }
            body { line-height: 1.5; color: #000000; font-size: 14px; }
            /* The Cytoscape Web container must have its dimensions set. */
            #cytoscapeweb { width: 100%; height: 80%; }
            #note { width: 100%; text-align: center; padding-top: 1em; }
            .link { text-decoration: underline; color: #0B94B1; cursor: pointer; }
        </style>
    </head>
    
    <body>
        <div id="cytoscapeweb">
            Cytoscape Web will replace the contents of this div with your graph.
        </div>
        <div id="note">
            Click node or edge for details
        </div>
    </body>
    
</html>
"""

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

class rgb:
    def __init__(self,r,g,b):
        self.r = int(round(r))
        self.g = int(round(g))
        self.b = int(round(b))
        
        if self.r > 255:
            self.r = 255
        elif self.r < 0:
            self.r = 0
        if self.g > 255:
            self.g = 255
        elif self.g < 0:
            self.g = 0
        if self.b > 255:
            self.b = 255
        elif self.b < 0:
            self.b = 0
        
    def tohex(self):
        r = self.r
        g = self.g
        b = self.b
        hexchars = "0123456789ABCDEF"
        return "#" + hexchars[r/16]+hexchars[r%16]+hexchars[g/16]+hexchars[g%16]+hexchars[b/16]+hexchars[b%16]

def getColor(val, minVal = -8, maxVal = 8):
    minColor = rgb(0, 0, 255)
    maxColor = rgb(255, 0, 0)
    zeroColor = rgb(255, 255, 255)
    try:
        fval = float(val)
    except ValueError:
        col = rgb(200, 200, 200)
        return col.tohex()
    if fval < 0.0:
        if fval < minVal:
            fval = 1.0
        else:
            fval = fval/minVal
        col = minColor
    else:
        if fval > maxVal:
            fval = 1.0
        else:
            fval = fval/maxVal
        col = maxColor
    r = fval*float(col.r-zeroColor.r)+zeroColor.r
    g = fval*float(col.g-zeroColor.g)+zeroColor.g
    b = fval*float(col.b-zeroColor.b)+zeroColor.b
    col = rgb(r,g,b)
    
    return(col.tohex())

def getColorbyType(feature, typeData):
    if typeData[feature] == "Basal":
        col = rgb(255, 0, 0)
    elif typeData[feature] == "Classical":
        col = rgb(0, 0, 255)
    elif typeData[feature] == "Primitive":
        col =  rgb(0, 255, 0)
    elif typeData[feature] == "Secretory":
        col = rgb(200, 0, 200)
    else:
        col = rgb(255, 255, 255)
    return(col.tohex())

def getSize(val, minVal = -10, maxVal = 10):
    fval = float(val)
    if fval < 0.0:
        if fval < minVal:
            size = 100
        else:
            size = 100*(fval/minVal)
    else:
        if fval > maxVal:
            size = 100
        else:
            size = 100*(fval/maxVal)
    return(size)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "n:sq")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        print "incorrect number of arguments"
        usage(1)
    
    feature = args[0]
    destDir = args[1]
    
    mutTable = False
    scoreTable = False
    customImage = False
    global verbose, noteText
    for o, a in opts:
        if o == "-n":
            noteText = a
        elif o == "-s":
            scoreTable = True
        elif o == "-q":
            verbose = False
    if destDir.endswith("/"):
        destDir = destDir.rstrip("/")
    
    ## check structure
    assert os.path.exists("LABEL.NA")
    assert os.path.exists("TYPE.NA")
    assert os.path.exists("%s_SCORE.NA" % (feature))
    assert os.path.exists("%s" % (feature))
    if feature.startswith("disc-plot"):
        mutTable = True
        mutGene = re.sub("disc-plot", "", feature)
    if os.path.exists("%s/img" % (feature)):
        customImage = True

    ## identify nets with feature
    sifFile = None
    for i in os.listdir("%s" % (feature+"/")):
        if i.endswith(netExtension):
            sifFile = feature+"/"+i
            break
    (allNodes, allInteractions) = mPathway.rSIF(sifFile)
    idMap = dict()
    nodeMap = dict()
    for i, node in enumerate(allNodes.keys()):
        idMap[node] = i+1
        nodeMap[i+1] = node
    labelMap = mData.r2Col("LABEL.NA", delim = " = ", header = True)
    typeMap = mData.r2Col("TYPE.NA", delim = " = ", header = True)
    scoreMap = mData.r2Col("%s_SCORE.NA" % (feature), delim = " = ", header = True)
    nodes = nodeMap.keys()
    nodes.sort()
    
    ## create graphml structure
    graphmlContent = """<graphml>\\
                    <key id="name" for="node" attr.name="name" attr.type="string"/>\\
                    <key id="label" for="node" attr.name="label" attr.type="string"/>\\
                    <key id="type" for="node" attr.name="type" attr.type="string"/>\\
                    <key id="color" for="node" attr.name="color" attr.type="string"/>\\
                    <key id="size" for="node" attr.name="size" attr.type="double"/>\\
                    <key id="score" for="node" attr.name="score" attr.type="double"/>\\
                    <key id="image" for="node" attr.name="image" attr.type="string"/>\\
                    <key id="interaction" for="edge" attr.name="interaction" attr.type="string"/>\\
                    <graph edgedefault="directed">\\
                    """
    nodeVals = []
    for i in nodes:
        if nodeMap[i] == "__DISCONNECTED__":
            nodeName = re.sub("'", "", nodeMap[i])
            nodeLabel = re.sub("'", "", nodeMap[i])
            nodeType = "protein"
            nodeColor = "#FFFFFF"
            nodeSize = 25
            nodeScore = 0
            nodeImage = ""
        else:
            try:
                nodeVals.append(abs(float(scoreMap[nodeMap[i]])))
            except ValueError:
                pass
            nodeName = re.sub("'", "", nodeMap[i])
            nodeLabel = re.sub("'", "", labelMap[nodeMap[i]])
            if nodeLabel == "NA":
                nodeLabel == ""
            nodeType = typeMap[nodeMap[i]]
            nodeColor = getColor(scoreMap[nodeMap[i]])
            nodeSize = getSize(scoreMap[nodeMap[i]])
            nodeScore = scoreMap[nodeMap[i]]
            nodeImage = ""
            # & (typeMap[nodeMap[i]] == "protein")
            if (customImage):
                if os.path.exists("%s/img/%s.png" % (feature, re.sub("[:/]", "_", nodeMap[i]))):
                    nodeImage = "img_%s/%s.png" % (feature, re.sub("[:/]", "_", nodeMap[i]))
                    nodeType = "plot"
                    nodeColor = nodeColor
                else:
                    nodeImage = ""
                    nodeType = nodeType
                    nodeColor = "#FFFFFF"
                    
        graphmlContent += """       <node id="%s">\\
                            <data key="name">%s</data>\\
                            <data key="label">%s</data>\\
                            <data key="type">%s</data>\\
                            <data key="color">%s</data>\\
                            <data key="size">%s</data>\\
                            <data key="score">%s</data>\\
                            <data key="image">%s</data>\\
                        </node>\\
                        """ % (i, nodeName, nodeLabel, nodeType, nodeColor, nodeSize, nodeScore, nodeImage)
    for i in allInteractions.keys():
        for j in allInteractions[i].keys():
            graphmlContent += """   <edge source="%s" target="%s">\\
                            <data key="interaction">%s</data>\\
                        </edge>\\
                        """ % (idMap[i], idMap[j], allInteractions[i][j])
    graphmlContent += """</graph>\\
                </graphml>\\
                """
    
    ## launch cytoscape
    if (not os.path.exists(destDir)):
        os.system("mkdir %s" % (destDir))
    f = open("%s/%s.html" % (destDir, feature), "w")
    f.write(htmlHead+graphmlContent+htmlTail)
    f.close()
    if customImage:
        os.system("cp -r %s/img %s/img_%s" % (feature, destDir, feature))
    
    ## mut table
    if mutTable:
        os.system("cp *.pdf %s/img_%s" % (destDir, feature))
        if os.path.exists("%s/stats.tab" % (destDir)):
            f = open("%s/stats.tab" % (destDir), "a")
        else:
            f = open("%s/stats.tab" % (destDir), "w")
            f.write("Feature\tCircle\tNonMut\tMut\tStrengthScore\tMutSig\tAUC\tSignal\tBackground\tSampleTable\n")
        s = open("sig.stats", "r")
        stats = re.split("\t", s.readline().rstrip("\t\r\n"))
        s.close()
        nnonmut = stats[1]
        nmut = stats[2]
        testAUC = stats[3]
        strengthScore = stats[4]
        pMutSig = stats[5]
        circleLink = htmlLink % ("%s/img_%s/%s.png" % (re.split("/", destDir)[-1], feature, mutGene), "png")
        aucLink = htmlLink % ("%s/img_%s/%s.auc.pdf" % (re.split("/", destDir)[-1], feature, mutGene), "pdf (%s)" % (testAUC))
        signalLink = htmlLink % ("%s/img_%s/%s.signal.pdf" % (re.split("/", destDir)[-1], feature, mutGene), "pdf")
        backgroundLink = htmlLink % ("%s/img_%s/%s.background.pdf" % (re.split("/", destDir)[-1], feature, mutGene), "pdf")
        tableLink = htmlLink % ("%s/table-%s.html" % (re.split("/", destDir)[-1], mutGene), "table")
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (feature, circleLink, nnonmut, nmut, strengthScore, pMutSig, aucLink, signalLink, backgroundLink, tableLink))
        f.close()
    
    ## score table
    elif scoreTable:
        if os.path.exists("%s/stats.tab" % (destDir)):
            f = open("%s/stats.tab" % (destDir), "a")
        else:
            f = open("%s/stats.tab" % (destDir), "w")
            f.write("id\tMAX\tUQS\ttot_nodes\tnote\n")
        f.write("%s\t%s\t%s\t%s\t%s\n" % (feature, max(nodeVals), mCalculate.quartiles(nodeVals)[2], len(allNodes.keys()), noteText))
        f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
