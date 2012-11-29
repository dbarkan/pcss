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
logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestPsipred(pcssTests.TestSequenceFeatures):

    def setUp(self):
        self.globalSetup("testConfig/testPcssConfig.txt")
        self.fileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.runner.pdh)
        self.sequenceFeatureReader = pcssFeatureHandlers.PsipredReader(self.fileHandler)
        self.sequenceFeatureRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)
        self.seqData = pcssTests.PsipredData()

    def getSeqFeatureCall(self, position):
        self.proteins[0].psipredProteinCalls.getSequenceFeatureCall(position)

    def getSeqFeatureCallMethod(self):
        return self.proteins[0].psipredProteinCalls.getSequenceFeatureCall
        

    def setBadCommandData(self):
        self.fileHandler.rootDataDir = os.path.join(self.pcssConfig["home_test_directory"], 
                                                            "testInput/psipredErrors/badCommand/psipredResults/")

        self.fileHandler.sequenceCmd = os.path.join(self.pcssConfig["home_test_directory"], 
                                                            "testInput/psipredErrors/badCommand/runPsipred/runpsipred")

    def getCallString(self):
        return self.proteins[0].psipredProteinCalls.makeFullCallString()

    def getProcessResultFileMethod(self):
        return self.proteins[0].processPsipred

    def processResultFile(self):
        self.proteins[0].processPsipred(self.sequenceFeatureReader, self.sequenceFeatureRunner)

if __name__ == '__main__':
    unittest.main()
