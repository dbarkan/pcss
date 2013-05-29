import unittest
import configobj
import pcssTools
import pcssPeptide
import os
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import logging
import sys
import pcssTests

class DisopredData:
    def __init__(self):

        self.stringFeatureValue = "OOOOOOOO"
        self.scoreFeatureValue = "0.014, 0.007, 0.007, 0.005, 0.007, 0.003, 0.003, 0.003"
                
    def getExpectedFullStringResult(self):
        return "DDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDDDDDDD"

class TestDisopred(pcssTests.TestSequenceFeatures):

    def setupSpecificTest(self):
        self.readProteins()
        self.seqData = DisopredData()
        self.fileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.runner.pdh)
        self.sequenceFeatureReader = pcssFeatureHandlers.DisopredReader(self.fileHandler)
        self.sequenceFeatureRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)
        self.errorDirName = "disopredErrors"
        self.name = "disopred"

    def getSeqFeatureCallMethod(self):
        return self.proteins[0].disorderProteinCalls.getSequenceFeatureCall

    def setBadCommandData(self):
        self.internalConfig["root_disopred_dir"] = os.path.join(self.getErrorInputFile("badCommand"), "disopredResults")
        self.fileHandler.sequenceCmd = os.path.join(self.getErrorInputFile("badCommand"), "runDisopred", "rundisopred")

    def setupLongRootDir(self):
        self.internalConfig["root_disopred_dir"] = "%(pcss_directory)s/test/testFileOutput/disopredResults/"

    def getCallString(self):
        return self.proteins[0].disopredProteinCalls.makeFullCallString()

    def getProcessResultFileMethod(self):
        return self.proteins[0].processDisopred

    def processResultFile(self):
        self.proteins[0].processDisopred(self.sequenceFeatureReader, self.sequenceFeatureRunner)

if __name__ == '__main__':
    unittest.main()
