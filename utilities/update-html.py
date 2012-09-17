#!/usr/bin/env python
"""update-html.py: 

Usage:
  update-html.py [options]

Options:
  -q            run quietly
"""
import os, os.path, sys, getopt, re
import mData

verbose = True

htmlHead = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html>

    <head>
        <title>CytoscapeWeb - UCSC PARADIGM Subnets</title>
        <script type="text/javascript" src="js/jquery-1.6.1.min.js"></script>
        <script type="text/javascript" src="js/jquery.tablesorter.js"></script>
        <script type="text/javascript">
            $(document).ready(function() 
                { 
                    $("#subnetTable").tablesorter(); 
                }
            );
        </script>
    </head>
    <body>
"""

htmlIndexItem = """        <a href="%s">%s</a><br>
"""

htmlCategory = """    <h1>%s</h1>
"""

htmlTail = """    </body>
</html>
"""

def usage(code = 0):
    print __doc__
    if code != None:
        sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "i:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 0:
        print "incorrect number of arguments"
        usage(1)
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## expand
    currentFiles = os.listdir(".")
    currentFiles.sort()
    d = open("index.html", "w")
    d.write(htmlHead)
    
    ## remove old files
    for i in currentFiles:
        if i.endswith(".html"):
            if i != "index.html":
                os.system("rm -f %s" % (i))
    
    ## work on new files
    for i in currentFiles:
        if i.endswith(".html"):
            continue
        elif i.endswith(".py"):
            continue
        elif i.endswith("js"):
            continue
        elif i.endswith("swf"):
            continue
        else:
            log("Working %s ...\n" % (i))
            d.write(htmlIndexItem % ("%s/stats.html" % (i), i))
    d.write(htmlTail)
    d.close()

if __name__ == "__main__":
    main(sys.argv[1:])