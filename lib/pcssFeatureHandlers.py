import sys
import os
from Bio import SeqIO
import itertools
import pcssTools
import logging
import subprocess
import pcssErrors
import pcssFeatures
import tempfile
import shutil
import time
log = logging.getLogger("pcssFeatureHandlers")


class SequenceFeatureFileHandler:
    
    """Base class for handling sequence-feature specific file and directory information"""

    def getSequenceFeatureFile(self, modbaseSeqId):
        """Return the full path file name for the sequence feature containing results for input sequence ID"""
        return os.path.join(self.getSequenceFeatureDir(modbaseSeqId),
                            "%s.%s" % (modbaseSeqId, self.outputFileSuffix))

    def getSequenceFeatureDir(self, modbaseSeqId):
        """Return the full path directory where input sequence id results will be written"""
        return os.path.join(self.getRootDataDir(), self.pdh.getTwoLetterOutputDir(modbaseSeqId))

    def outputFileExists(self, modbaseSeqId):
        return os.path.exists(self.getSequenceFeatureFile(modbaseSeqId))

class DisopredFileHandler(SequenceFeatureFileHandler):
    """Class providing Disopred-specific data for running algorithm and processing results"""
    def __init__(self, pcssConfig, pdh):
        self.outputFileSuffix = "diso"
        self.sequenceCmd = pcssConfig["run_disopred_command"]
        self.pdh = pdh
        self.pcssConfig = pcssConfig
        self.name = "disopred"

    def getRootDataDir(self):
        return self.pcssConfig["root_disopred_dir"]

    def getName(self):
        return self.name

    def getCommandException(self, msg):
        return pcssErrors.PcssGlobalException(msg)
    
class PsipredFileHandler(SequenceFeatureFileHandler):
    """Class providing Psipred-specific data for running algorithm and processing results"""
    def __init__(self, pcssConfig, pdh):
        self.outputFileSuffix = "ss2"
        self.sequenceCmd = pcssConfig["run_psipred_command"]
        self.pdh = pdh
        self.pcssConfig = pcssConfig
        self.name = "psipred"

    def getRootDataDir(self):
        return self.pcssConfig["root_psipred_dir"]
                
    def getName(self):
        return self.name

    def getCommandException(self, msg):
        return pcssErrors.PcssGlobalException(msg)

class SequenceFeatureRunner:
    
    """Class that runs algorithm to calculate a sequence feature for a pcssProtein"""

    def __init__(self, sequenceFileHandler):
        self.sfh = sequenceFileHandler

    def runSequenceFeature(self, pcssProtein):
        """Call command for running shell script to run the sequence feature algorithm. 

        Script runs in a temporary directory and copies results to permanent output location, defined by parameters
        """
        tempDirectory = pcssProtein.pcssRunner.pdh.makeTempDirectory()

        try:
            self.runSequenceSubprocess(pcssProtein, tempDirectory.tempDir)
        except pcssErrors.ProteinException as e:
            tempDirectory.changeBack()
            raise e
        except pcssErrors.PcssGlobalException as e:
            tempDirectory.changeBack()
            raise e
        outputFile = self.getSequenceFeatureOutputFile(pcssProtein)
        destinationDir = self.makeSeqFeatureDirectory(pcssProtein)
        
        pcssProtein.pcssRunner.pdh.moveFile(tempDirectory.tempDir, outputFile, destinationDir)
        tempDirectory.changeBack()
        

    def getSequenceFeatureOutputFile(self, pcssProtein):
        return "%s.%s" % (pcssProtein.modbaseSequenceId, self.sfh.outputFileSuffix)

    def runSequenceSubprocess(self, pcssProtein, cwd):
        """Run system command that runs disopred / psipred algorithm"""
        fastaFile = pcssProtein.writeSequenceToFasta(cwd)
        print "running command %s" % self.sfh.sequenceCmd
        commandList = self.sfh.sequenceCmd.split(" ")
        commandList.append(fastaFile)
        
        output = self.sfh.pdh.runSubprocess(commandList, checkStdError=False)
        
        outputFile = self.getSequenceFeatureOutputFile(pcssProtein)
        if (not os.path.exists(outputFile)):

            raise self.sfh.getCommandException("Finished %s command\nOutput file not where expected (looked for %s).\n "
                                        "Note that disopred/psipred can take hours to run on long sequences; "
                                        "Your job might have timed out. Command output: \n%s" %
                                        (self.sfh.sequenceCmd, outputFile, output[1]))
    
    def makeSeqFeatureDirectory(self, pcssProtein):
        """Make directory where sequence feature result will be stored (two letter prefix)"""
        destinationDir = self.sfh.getSequenceFeatureDir(pcssProtein.modbaseSequenceId)

        self.sfh.pdh.runSubprocess(["mkdir", "-p", destinationDir])
        self.sfh.pdh.sleepUntilDone(destinationDir, predicate=self.sfh.pdh.fileDoesNotExist)

        return destinationDir

class SequenceFeatureReader:
    
    """Abstract class for reading Sequence Feature Result; needs to be sublcassed"""

    def readResult(self, modbaseSeqId):
        print "Error: sequenceFeatureReader should only be implemented as a subclass"
        sys.exit()
        
class DisopredReader(SequenceFeatureReader):

    """Class for reading Disopred result file and creating a SeuqenceFeatureCalls to store it"""

    def __init__(self, disopredFileHandler):
        self.sfh = disopredFileHandler

    def readResult(self, modbaseSeqId):
        """Read .diso result file and store all data as a SequenceFeatureCalls"""
        sequenceFile = self.sfh.getSequenceFeatureFile(modbaseSeqId)
        reader = pcssTools.PcssFileReader(sequenceFile)
        lines = reader.getLines()
        dpc = pcssFeatures.DisopredSequenceFeatureCallSet()
        for i in range(4, len(lines)): #first four lines aren't commented but need to be skipped
            disorderCall = pcssFeatures.DisorderResidueCall(lines[i])
            dpc.addSequenceFeatureCall(disorderCall)
        
        return dpc

class PsipredReader(SequenceFeatureReader):

    """Class for reading Psipred result file and creating a SequenceFeatureCalls to store it"""

    def __init__(self, psipredFileHandler):
        self.sfh = psipredFileHandler

    def readResult(self, modbaseSeqId):
        """Read .ss2 result file and store all data as a SequenceFeatureCalls"""
        sequenceFile = self.sfh.getSequenceFeatureFile(modbaseSeqId)
        reader = pcssTools.PcssFileReader(sequenceFile)
        lines = reader.getLines()
        ppc = pcssFeatures.PsipredSequenceFeatureCallSet()
        for i in range(len(lines)): 
            psipredCall = pcssFeatures.PsipredResidueCall(lines[i])
            ppc.addSequenceFeatureCall(psipredCall)
        return ppc
    
