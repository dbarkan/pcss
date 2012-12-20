import sys
import re

import os
import pcssTools

if (len(sys.argv) < 3):
    print "Usage: python cleanModelTable.py <inputFile> <outputFile>"
    print "Remove last empty tab character from <inputFile> for each line and write to <outputFile>"
    sys.exit()

inputFile = sys.argv[1]
outputFile = sys.argv[2]
reader = pcssTools.PcssFileReader(inputFile)
outputFh = open(outputFile, 'w')
lines = reader.getLines()
for line in lines:
    cols = line.split('\t')
    if (cols[-1] != ""):
        print "error: got line with non-empty last column: %s" % line
        sys.exit()
    outputLine = '\t'.join(cols[0:len(cols) -1])
    outputFh.write("%s\n" % outputLine)

