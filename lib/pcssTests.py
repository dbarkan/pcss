import unittest
import configobj
import pcssIO
import pcssTools
import pcssPeptide
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import logging
import sys
import os
import pcssTests
logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestSequenceFeaturesLong(unittest.TestCase):
    def globalSetup(self, configFile):
        configSpecFile = "testConfig/testConfigSpec.txt"        

        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])

class TestSequenceFeatures(unittest.TestCase):

    def updateConfig(self, name, value):
        self.pcssConfig['name'] = value

    def globalSetup(self, configFile):

        configSpecFile = "testConfig/testConfigSpec.txt"        

        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])

    def test_parse_result_file(self):
        self.processResultFile()

        callString = self.getCallString()

        self.assertEquals(callString, self.seqData.getExpectedFullStringResult())
        
        peptide = self.proteins[0].peptides[17]

        self.assertEquals(peptide.attributes[self.seqData.stringFeatureName].getValueString(), self.seqData.stringFeatureValue)
        self.assertEquals(peptide.attributes[self.seqData.scoreFeatureName].getValueString(), self.seqData.scoreFeatureValue)

        self.proteins[0].setPeptide(pcssPeptide.PcssPeptide("FAKE", 1000, 1003, self.runner))
        self.processException(self.proteins[0].peptides[1000], self.seqData.stringFeatureName, self.seqData.peptideNotFoundCode, self.processResultFile)

    def processException(self, peptide, errorAttribute, exceptionCode, function, *args):
        function(*args)

        self.assertEquals(peptide.getAttributeOutputString(errorAttribute), exceptionCode)
        self.assertRegexpMatches(peptide.getAttributeOutputString("peptide_errors"), exceptionCode)
        
    def test_found_disopred_file(self):
        
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"        

        pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)

        runner = pcssTools.PcssRunner(pcssConfig)

        self.assertTrue(self.fileHandler.outputFileExists("76c3a409540532138c6b44bde9e4d248MDDRDENQ"))

    def test_residue_mismatch(self):
        
        self.fileHandler.rootDataDir = self.seqData.residueMismatchDir
        self.processException(self.proteins[0].peptides.values()[0], self.seqData.stringFeatureName, self.seqData.proteinMismatchCode, self.processResultFile)

    def test_unexpected_call(self):
        self.fileHandler.rootDataDir = self.seqData.badCallDir
        self.processException(self.proteins[0].peptides.values()[0], self.seqData.stringFeatureName, self.seqData.badCallCode, self.processResultFile)

    def test_command_error(self):
        initialCwd = os.getcwd()
        self.setBadCommandData()
        self.processException(self.proteins[0].peptides.values()[0], self.seqData.stringFeatureName, self.seqData.badCommandCode, self.processResultFile)
        self.assertEquals(os.getcwd(), initialCwd)

    def test_bad_line(self):
        self.fileHandler.rootDataDir = self.seqData.badLineDir
        self.processException(self.proteins[0].peptides.values()[0], self.seqData.stringFeatureName, self.seqData.badLineCode, self.processResultFile)

class DisopredData:
    def __init__(self):
        self.name = "disopred"

        self.stringFeatureName = "disopred_string_feature"
        self.stringFeatureValue = "OOOOOOOO"
        self.scoreFeatureName =  "disopred_score_feature"
        self.scoreFeatureValue = "0.016, 0.01, 0.007, 0.003, 0.008, 0.004, 0.003, 0.002"

        self.residueMismatchDir = "testInput/disopredErrors/residueMismatch/"
        self.badCallDir = "testInput/disopredErrors/badCall/"
        self.badLineDir = "testInput/disopredErrors/badLine/"

        self.peptideNotFoundCode = "peptide_disopred_peptide_not_found"
        self.proteinMismatchCode = "peptide_disopred_protein_mismatch"
        self.badCommandCode = "peptide_disopred_bad_command"
        self.badLineCode = "peptide_disopred_bad_line"
        self.badCallCode = "peptide_disopred_bad_call"

    def getExpectedFullStringResult(self):
        return "DDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDOOOOOOOOOOOOOOOOODDDDDDDDDDDDDDDDDDDDDDDD"

class PsipredData:
    def __init__(self):
        self.stringFeatureName = "psipred_string_feature"
        self.stringFeatureValue = "AAAAAAAA"
        self.scoreFeatureName =  "psipred_score_feature"
        self.scoreFeatureValue =  "0.256, 0.602, 0.049, 0.011, 0.008, 0.004, 0.003, 0.004"
        self.badLineDir = "testInput/psipredErrors/badLine/"
        self.residueMismatchDir = "testInput/psipredErrors/residueMismatch/"
        self.badCallDir = "testInput/psipredErrors/badCall/"
        self.peptideNotFoundCode = "peptide_psipred_peptide_not_found"
        self.proteinMismatchCode = "peptide_psipred_protein_mismatch"
        self.badCommandCode = "peptide_psipred_bad_command"
        self.badCallCode = "peptide_psipred_bad_call"
        self.badLineCode = "peptide_psipred_bad_line"

    def getExpectedFullStringResult(self):
        return "LLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALLLLLLAAAAAAAAAAAAAAALLAAAAAAAAAAAAAAAALLLLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALLLLLLLLLLBBBBBBLLLLBBBBLLLLLLLAAAAAAAAAAAAAAAAAAAAAAAALLLLLLLLLLAAAAAAAAAAAAALLAAAAAAAAAAAAAAAAAAALLLLLLAAAAAAAAAAAAAAAAAALLLLLLLLLLLLLLLLLLLLLLLLL"
