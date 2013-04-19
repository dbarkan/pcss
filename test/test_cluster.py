import sys
import os
import pcssTools
import pcssCluster
import pcssErrors
import unittest
import configobj

class TestCluster(unittest.TestCase):


    def setUp(self):
        configFile = "testConfig/svmTrainingBenchmarkConfig.txt"
        tempPcssConfig = configobj.ConfigObj(configFile)
        configSpecFile = tempPcssConfig["user_config_spec_file"]

        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")

        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        
    def dtest_preprocess(self):

        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")

        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        self.runner.setJobDirectory(os.path.join(self.runner.pcssConfig["run_directory"], "developClusterJob"))
        seqDivider = pcssCluster.SeqDivider(self.runner)
        seqDivider.divideSeqsFromFasta()
        
        seqDivider.csg.makeFullSgeScript() #needs update

    def test_missing_input_annotation_file(self):
        self.runner.internalConfig["annotation_output_file"] = "fake"
        cb = pcssCluster.ClusterBenchmarker(self.runner)
        self.assertRaises(pcssErrors.PcssGlobalException, cb.prepareTrainingBenchmarkRun)

if __name__ == '__main__':
    unittest.main()

