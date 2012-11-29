import unittest
import configobj
import pcssTools
import pcssPeptide


class TestReadInput(unittest.TestCase):
    
    def test_scan_peptides(self):
        configFile = "testConfig/testPcssConfig.txt"

        pcssConfig = configobj.ConfigObj(configFile)

        pcssRunner = pcssTools.PcssRunner(pcssConfig)

        spi = pcssPeptide.ScanPeptideImporter(pcssRunner)
        proteins = spi.readInputFile(pcssRunner.pcssConfig['fasta_file'])
        
        self.assertEqual(proteins[0].modbaseSequenceId, "39d0244382b7121b3d6645150e8b77feMRVTSLTA")
        self.assertEqual(proteins[0].uniprotId, "P30498")
        self.assertEqual(proteins[0].peptides.values()[0].getOutputString(), "Start: 49 end 56 sequence GYVDDTQF\n")
        self.assertEqual(len(proteins[0].peptides.values()), 4)
        

if __name__ == '__main__':
    unittest.main()
