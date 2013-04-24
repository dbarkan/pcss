import unittest
import os
import configobj
import pcssTools
import pcssPeptide
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import logging
import sys
import pcssTests


class PsipredData:
    def __init__(self):
        self.stringFeatureValue = "AAAAAAAA"
        self.scoreFeatureValue =  "0.454, 0.624, 0.122, 0.017, 0.01, 0.004, 0.002, 0.004"

    def getExpectedFullStringResult(self):
        return "LLLAAAAAAAAAAAAAAAAAAAAAAAAAAAALLLLLLLLAAAAAAAAAAAAAALLLAAAAAAAAAAAAAAALLLLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALLLLLLLLLLBBBAAAAAAAAAAAAAAALLLLAAAAAAAAAAAAAAAAAAAAAAALLLLLAAAAAAAAAAAAAAAAALLLAAAAAAAAAAAAAAAAAAAAALLLAAAAAAAAAAAAAAAAAAALLLLLLLLLLLLLLLLLLLLLLLLL"

class TestPsipred(pcssTests.TestSequenceFeatures):

    def setupSpecificTest(self):
        self.readProteins()
        self.seqData = PsipredData()
        self.fileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.runner.pdh)
        self.sequenceFeatureReader = pcssFeatureHandlers.PsipredReader(self.fileHandler)
        self.sequenceFeatureRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)

        self.errorDirName = "psipredErrors"
        self.name = "psipred"

    def getSeqFeatureCall(self, position):
        self.proteins[0].psipredProteinCalls.getSequenceFeatureCall(position)

    def getSeqFeatureCallMethod(self):
        return self.proteins[0].psipredProteinCalls.getSequenceFeatureCall
        
    def setBadCommandData(self):
        self.pcssConfig["root_psipred_dir"] = os.path.join(self.getErrorInputFile("badCommand"), "psipredResults")
        self.fileHandler.sequenceCmd = os.path.join(self.getErrorInputFile("badCommand"), "runPsipred", "runpsipred")

    def getCallString(self):
        return self.proteins[0].psipredProteinCalls.makeFullCallString()

    def getProcessResultFileMethod(self):
        return self.proteins[0].processPsipred

    def processResultFile(self):
        self.proteins[0].processPsipred(self.sequenceFeatureReader, self.sequenceFeatureRunner)

    def setupLongRootDir(self):
        self.pcssConfig["root_psipred_dir"] = "%(pcss_directory)s/test/testFileOutput/psipredResults/"

        

if __name__ == '__main__':
    unittest.main()
