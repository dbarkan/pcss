import unittest
import configobj
import pcssTools

import pcssPeptide
import pcssTests
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import shutil
import logging
import sys
import time
import os
logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestSequenceFeaturesLong(pcssTests.TestSequenceFeaturesLong):

    def dtest_psipred_long(self):
        self.globalSetup("testConfig/testPcssConfigLong.txt")
        self.fileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.runner.pdh)
        psipredReader = pcssFeatureHandlers.PsipredReader(self.fileHandler)
        psipredRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)

        psipredData = pcssTests.PsipredData()
        self.proteins[0].processPsipred(psipredReader, psipredRunner)

        psipredCallString = self.proteins[0].psipredProteinCalls.makeFullCallString()
        self.assertEquals(psipredCallString, psipredData.getExpectedFullStringResult())

        peptide = self.proteins[0].peptides[17]
        self.assertEquals(peptide.features["psipred_string_feature"].getValueString(), psipredData.stringFeatureValue)
        self.assertEquals(peptide.features["psipred_score_feature"].getValueString(), psipredData.scoreFeatureValue)
        self.assertRaises(pcssErrors.FeatureException, self.proteins[0].psipredProteinCalls.getSequenceFeatureCall, 1000)
        time.sleep(5)
        shutil.rmtree(self.fileHandler.getSequenceFeatureDir(self.proteins[0].modbaseSequenceId))

    def test_disopred(self):
        print  time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())
        self.globalSetup("testConfig/testPcssConfigLong.txt")
        self.fileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.runner.pdh)
        disopredReader = pcssFeatureHandlers.DisopredReader(self.fileHandler)
        disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(self.fileHandler)

        disopredData = pcssTests.DisopredData()
        self.proteins[0].processDisopred(disopredReader, disopredRunner)

        disopredCallString = self.proteins[0].disorderProteinCalls.makeFullCallString()
        self.assertEquals(disopredCallString, disopredData.getExpectedFullStringResult())

        peptide = self.proteins[0].peptides[17]
        self.assertEquals(peptide.features["disorder_string_feature"].getValueString(), disopredData.stringFeatureValue)
        self.assertEquals(peptide.features["disorder_score_feature"].getValueString(), disopredData.scoreFeatureValue)
        self.assertRaises(pcssErrors.FeatureException, self.proteins[0].disorderProteinCalls.getSequenceFeatureCall, 1000)
        time.sleep(5)
        shutil.rmtree(self.fileHandler.getSequenceFeatureDir(self.proteins[0].modbaseSequenceId))
        print  time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())
    
if __name__ == '__main__':
    unittest.main()