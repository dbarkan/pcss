import sys
import os
import pcssTools
import pcssCluster
import unittest
import configobj

class TestCluster(unittest.TestCase):

    def test_preprocess(self):

        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")

        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        seqDivider = pcssCluster.SeqDivider(self.runner)
        seqDivider.divideSeqsFromFasta(self.pcssConfig['fasta_file'])
        
        seqDivider.makeFullSgeScript()

if __name__ == '__main__':
    unittest.main()
