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

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestDisopred(pcssTests.TestSequenceFeatures):

    def setUp(self):
        self.globalSetup("testConfig/testPcssConfig.txt")
        self.seqData = pcssTests.DisopredData()
        self.fileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.runner.pdh)
        self.sequenceFeatureReader = pcssFeatureHandlers.DisopredReader(self.fileHandler)
        self.sequenceFeatureRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)


    def getSeqFeatureCallMethod(self):
        return self.proteins[0].disorderProteinCalls.getSequenceFeatureCall

    def setBadCommandData(self):
        self.fileHandler.rootDataDir = os.path.join(self.pcssConfig["home_test_directory"], 
                                                            "testInput/disopredErrors/badCommand/disopredResults/")

        self.fileHandler.sequenceCmd = os.path.join(self.pcssConfig["home_test_directory"], 
                                                            "testInput/disopredErrors/badCommand/runDisopred/rundisopred")

    def getCallString(self):
        return self.proteins[0].disopredProteinCalls.makeFullCallString()

    def getProcessResultFileMethod(self):
        return self.proteins[0].processDisopred

    def processResultFile(self):
        print "running process result file"
        self.proteins[0].processDisopred(self.sequenceFeatureReader, self.sequenceFeatureRunner)
        

if __name__ == '__main__':
    unittest.main()
