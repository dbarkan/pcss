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
import pcssSvm
import pcssIO
from validate import Validator

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestRunner(unittest.TestCase):

    def setUp(self):
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)

    def clearErrorFiles(self):
        self.runner.pdh.runSubprocess(["rm", self.runner.pdh.getFullOutputFile("errors.out")], False)

    def test_pcss_error(self):
        self.runner = pcssTools.SvmApplicationRunner(self.pcssConfig)
        self.pcssConfig["peptide_importer_type"] = "defined"
        self.pcssConfig['fasta_file'] = "testInput/ioErrors/peptideMismatchFasta.txt"
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getFullOutputFile("errors.out")))

    def test_svm_runner(self):
        
        self.runner = pcssTools.SvmApplicationRunner(self.pcssConfig)

        self.runner.execute()

if __name__ == '__main__':
    unittest.main()
