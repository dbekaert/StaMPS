#!/usr/bin/python

# Code for extracting heading angle from TSX leader file
#
# Adapted from tsx_dump_header2doris.py
#
# Karsten Spaans, July 2010
# ===========================================================================
# ===========================================================================

from lxml import etree
import string, time, sys
#import xml.etree.ElementTree as ElementTree
#import types

def usage():
    print '\nUsage: python tsx_dump_header2doris.py tsx_XML_product > outputfile'
    print '  where tsx_XML_product is the input filename'
#    print '        outputfile      is the output DORIS resultfile'

try:
    inputFileName  = sys.argv[1]
#    outputFileName = sys.argv[2]
#    outStream      = open(outputFileName,'w')
except:
    print 'Unrecognized input'
    usage()
    sys.exit(1)

# input file
#inputFileName  = '04091/04091.xml'
#
#outputFileName = 'inputFileName'


inTree = etree.parse(inputFileName)

queryList = {\
             'heading'  : './/sceneInfo/headingAngle',\
	    }

# temp variables and parameters
container     = {}
# containerTemp = {}        
events        = ('end',)

# functions : not sure if needed
def fast_iter_string(context):
    for event,elem in context:
        return elem.text

# works with lists
def fast_iter_list(context,tag=''):
    for event,elem in context:
        return elem.iterchildren(tag=tag).next().text

def hms2sec(hmsString,convertFlag='int'):
    # input hmsString syntax: XX:XX:XX.xxxxxx
    secString = int(hmsString[0:2])*3600 + \
        int(hmsString[3:5])*60 + \
        float(hmsString[6:])
    if convertFlag == 'int' :
        return int(secString)
    elif convertFlag == 'float' :
        return float(secString)
    else:
        return int(secString)


for key in queryList.keys():

    try:
        vars()[key];
    except KeyError or NameError:
        vars()[key] = [];

    queryListKey = queryList[key]
    for nodes in inTree.findall(queryListKey):
#    for nodes in inTree.findall(queryList[key]):
        vars()[key].append(nodes.text)

    container[key] = vars()[key]


# ---------------------------------------------------------------------------------------------------------

print('HEADING 					%s' % container['heading'][0])
