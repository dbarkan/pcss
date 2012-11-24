import sys
import os
import re
import subprocess

class PcssRunner:

    """Class for running PCSS algorithms'. Handles directory structure and some file conventions, among other tools"""

    def __init__(self, pcssConfig):
        """Create new PTP runner

        Construction does the following:
        1. Creates the 'run directory' which is where all output for this program goes. Run directory is created as $run_directory/$run_name
        where each variable specified in the parameter file. Thus if the user wants to run a program twice with different input, simply change
        the run_name parameter and nothing from previous runs will be overwritten
        """
        outputDir = pcssConfig["run_directory"]
        runName = pcssConfig["run_name"]
        print "output: %s run %s" % (outputDir, runName)
        fullOutputDir = os.path.join(outputDir, runName)
        mkdirProcess = subprocess.Popen(['mkdir', '-p', fullOutputDir], shell=False, stderr=subprocess.PIPE)
        self.outputDir =  os.path.join(outputDir, runName)

    def getFullOutputFile(self, fileName):
        """Easy way to get the full path of fileName where the path is the full run directory for this run"""
        return os.path.join(self.outputDir, fileName)

class PtFileReader:

    """Essentially python file object that provides commonly used processing functionality"""

    def __init__(self, fileName, skipBlanks=True, skipHash=True):
        """Read input file and save all lines after processing

        Similar to python File object, except we commonly will want to skip blank lines in the file
        (programs often break when getting an unexpected blank line) and also skip lines beginning with
        '#' as these are often comments as oposed to data.

        @param fileName full path file name to be processed
        @param skipBlanks set to true to ignore blank lines
        @param skipHash set to true to ignore lines beginning with '#'
        """
        fileHandle = open(fileName, 'r')
        lines = fileHandle.readlines()
        finalLines = []
        blankRe = re.compile('^\s*$')
        hashRe = re.compile('^\#')

        for line in lines:
            line = line.rstrip('\r\n')
            blankLine = blankRe.search(line)
            if (blankLine and skipBlanks):
                continue
            hashLine = hashRe.search(line)
            if (hashLine and skipHash):
                continue
            finalLines.append(line)
        self.finalLines = finalLines
        self.fileName = fileName

    def getLines(self):
        return self.finalLines
