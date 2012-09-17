#!/usr/bin/env python
"""pathmark-report.py: Generate cytoscapeweb sessions from a bundle (.sif and .NA files)

Usage:
  pathmark-report.py [options] bundleDir reportDir

Options:
  -t file[,file,...]         tab files containing statistics to include on the report
  -q                         run quietly
"""
## Written By: 
import os, os.path, sys, getopt, re
import mData, mPathway, mCalculate

verbose = True

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

htmlTableHead = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html>

    <head>
        <title>CytoscapeWeb - UCSC PARADIGM Subnets</title>
        <script type="text/javascript" src="../js/jquery-1.6.1.min.js"></script>
        <script type="text/javascript" src="../js/jquery.tablesorter.js"></script>
        <script type="text/javascript">
            $(document).ready(function() 
                { 
                    $("#htmlTable").tablesorter(); 
                }
            );
        </script>
    </head>
    <body>
"""

htmlTableHeader = """    <table id="htmlTable" class="tablesorter" border="1">
        <thead>
        <tr>
            <th>%s</th>
        </tr>
        </thead>
        <tbody>
"""

htmlTableItem = """        <tr>
            <td>%s</td>
        </tr>
"""

htmlTableTail = """        </tbody>
    </table>
    </body>
</html>
"""

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if (verbose):
        sys.stderr.write(msg)
    if (die):
        sys.exit(1)

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
        opts, args = getopt.getopt(args, "t:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        print "incorrect number of arguments"
        usage(1)
    
    featureDir = args[0].rstrip("/")
    featureName = featureDir.split("/")[-1]
    htmlDir = args[1].rstrip("/")
    
    tableFiles = []
    global verbose
    for o, a in opts:
        if o == "-t":
            tableFiles = a.split(",")
        elif o == "-q":
            verbose = False
    
    ## check structure
    assert os.path.exists("%s" % (featureDir))
    assert os.path.exists("%s/LABEL.NA" % (featureDir))
    assert os.path.exists("%s/TYPE.NA" % (featureDir))
    assert os.path.exists("%s/SCORE.NA" % (featureDir))
    
    ## read in pathway
    sifFile = None
    for i in os.listdir(featureDir):
        if i.endswith(".sif"):
            sifFile = "%s/%s" % (featureDir, i)
            break
    (nodes, interactions) = mPathway.rSIF(sifFile)
    
    ## construct mapping for integer id <-> string id
    idMap = {}
    nodeMap = {}
    for i, node in enumerate(nodes.keys()):
        idMap[node] = i+1
        nodeMap[i+1] = node
    labelMap = mData.r2Col("%s/LABEL.NA" % (featureDir), delim = " = ", header = True)
    typeMap = mData.r2Col("%s/TYPE.NA" % (featureDir), delim = " = ", header = True)
    scoreMap = mData.r2Col("%s/SCORE.NA" % (featureDir), delim = " = ", header = True)
    nodes = nodeMap.keys()
    nodes.sort()
    
    ## construct graphml
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
    
    ## add each node to graphml
    nodeVals = []
    for node in nodes:
        try:
            nodeVals.append(abs(float(scoreMap[nodeMap[node]])))
        except ValueError:
            pass
        except KeyError:
            pass
    if len(nodeVals) == 0:
        nodeVals = [-8, 8]
    for node in nodes:
        if nodeMap[node] == "__DISCONNECTED__":
            nodeName = re.sub("'", "", nodeMap[node])
            nodeLabel = re.sub("'", "", nodeMap[node])
            nodeType = "protein"
            nodeColor = "#FFFFFF"
            nodeSize = 25
            nodeScore = 0
            nodeImage = ""
        else:
            nodeName = re.sub("'", "", nodeMap[node])
            nodeLabel = re.sub("'", "", labelMap[nodeMap[node]])
            if nodeLabel == "NA":
                nodeLabel == ""
            nodeType = typeMap[nodeMap[node]]
            nodeColor = getColor(scoreMap[nodeMap[node]], minVal = min(nodeVals), 
                                 maxVal = max(nodeVals))
            nodeSize = getSize(scoreMap[nodeMap[node]], minVal = min(nodeVals),
                               maxVal = max(nodeVals))
            nodeScore = scoreMap[nodeMap[node]]
            nodeImage = ""
            if os.path.exists("%s/img" % (featureDir)):
                if os.path.exists("%s/img/%s.png" % (featureDir, re.sub("[:/]", "_", nodeMap[node]))):
                    nodeType = "plot"
                    nodeImage = "img_%s/%s.png" % (featureName, re.sub("[:/]", "_", nodeMap[node]))
                    nodeColor = nodeColor
                else:
                    nodeType = nodeType
                    nodeImage = ""
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
                        """ % (node, nodeName, nodeLabel, nodeType, nodeColor, nodeSize, nodeScore, nodeImage)
    
    ## add each connection to graphml
    for source in interactions.keys():
        for target in interactions[source].keys():
            graphmlContent += """   <edge source="%s" target="%s">\\
                            <data key="interaction">%s</data>\\
                        </edge>\\
                        """ % (idMap[source], idMap[target], interactions[source][target])
    graphmlContent += """</graph>\\
                </graphml>\\
                """
    
    ## create cytoscape-web html
    if (not os.path.exists(htmlDir)):
        assert(os.path.exists("/".join(htmlDir.split("/")[:-1])))
        os.system("mkdir %s" % (htmlDir))
    o = open("%s/%s.html" % (htmlDir, featureName), "w")
    o.write(htmlHead+graphmlContent+htmlTail)
    o.close()
    if os.path.exists("%s/img" % (featureDir)):
        os.system("cp -r %s/img %s/img_%s" % (featureDir, htmlDir, featureName))
        
    ## update table file
    tabHeader = ["cytoscapeweb"]
    tabLine = [htmlLink % ("%s.html" % (featureName), featureName)]
    for file in tableFiles:
        if file.endswith(".pdf"):
            (thisHeader, thisLine) = file.split(":")
            tabHeader.append(thisHeader)
            tabLine.append(htmlLink % ("pdf_%s/%s" % (featureName, thisLine.split("/")[-1]), "pdf"))
            if not os.path.exists("%s/pdf_%s" % (htmlDir, featureName)):
                os.system("mkdir %s/pdf_%s" % (htmlDir, featureName))
            os.system("cp %s %s/pdf_%s/" % (thisLine, htmlDir, featureName))
        else:
            f = open(file, "r")
            line = f.readline().rstrip()
            if line.startswith(">"):
                (thisHeader, thisLine) = line.split("\t")[1].split(":")
                tabHeader.append(thisHeader)
                tabLine.append(htmlLink % ("%s_%s.html" % (featureName, thisHeader), thisLine))
                o = open("%s/%s_%s.html" % (htmlDir, featureName, thisHeader), "w")
                o.write(htmlTableHead)
                o.write(htmlTableHeader % ("</th>\n            <th>".join(f.readline().lstrip("# ").rstrip().split("\t"))))
                for line in f:
                    o.write(htmlTableItem % ("</td>\n            <td>".join(line.rstrip().split("\t"))))
                o.write(htmlTableTail)
                o.close()
            elif line.startswith("#"):
                thisHeader = line.split("\t")
                thisLine = f.readline().rstrip().split("\t")
                assert(thisLine[0] == featureName)
                tabHeader += thisHeader[1:]
                tabLine += thisLine[1:]
            else:
                log("ERROR: Unrecognized file format for table addition: %s" % (file), die = True)
            f.close()
            
    if os.path.exists("%s/stats.tab" % (htmlDir)):
        o = open("%s/stats.tab" % (htmlDir), "a")
    else:
        o = open("%s/stats.tab" % (htmlDir), "w")     
        o.write("%s\n" % ("\t".join(tabHeader)))
    o.write("%s\n" % ("\t".join(tabLine)))
    o.close()
    
    ## rebuild html table
    f = open("%s/stats.tab" % (htmlDir), "r")
    o = open("%s/stats.html" % (htmlDir), "w")
    o.write(htmlTableHead)
    o.write(htmlTableHeader % ("</th>\n            <th>".join(f.readline().lstrip("# ").rstrip().split("\t"))))
    for line in f:
        o.write(htmlTableItem % ("</td>\n            <td>".join(line.rstrip().split("\t"))))
    o.write(htmlTableTail)
    o.close()
    f.close()
    ### pull this from update.html
    
if __name__ == "__main__":
    main(sys.argv[1:])
