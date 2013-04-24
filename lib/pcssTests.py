import unittest
import time
import configobj
import pcssIO
import pcssTools
import shutil
import pcssPeptide
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import logging
import sys
import os
import pcssTests

class PcssTest(unittest.TestCase):
    def setUp(self):
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"        
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.setupSpecificTest()

    def getProtein(self, proteinId, proteins):
        for protein in proteins:
            if (protein.modbaseSequenceId == proteinId):
                return protein

    def setSvmApplicationFileAttributes(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "svmApplicationFileAttributes.txt")

    def setAnnotationFileAttributes(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "annotationFileAttributes.txt")

    def setTrainingFileAttributes(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")

    def compareToExpectedOutput(self, observedFileName, expectedFileName, compareAlmostEqual=False):
        fullExpectedFileName = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "expectedOutput", 
                                            self.testName, "%s_expectedOutput.txt" % expectedFileName)
        print "compare %s and %s" % (observedFileName, fullExpectedFileName)
        if (compareAlmostEqual):
            self.compareFilesAlmostEqual(observedFileName, fullExpectedFileName)
        else:
            self.compareFiles(observedFileName, fullExpectedFileName, True)

    def setupSpecificTest(self):
        return

    def getErrorInputFile(self, fileName):
        return os.path.join(self.pcssConfig["home_test_directory"], "testInput", self.errorDirName, fileName)

    def compareFiles(self, firstFile, secondFile, sortLines=False):
        [firstLines, secondLines] = self.getLinesToCompare(firstFile, secondFile, sortLines)
        for (i, firstLine) in enumerate(firstLines):

            secondLine = secondLines[i]
            self.assertEquals(firstLine, secondLine)

    def compareFilesAlmostEqual(self, firstFile, secondFile, sortLines=False):
        [firstLines, secondLines] = self.getLinesToCompare(firstFile, secondFile, sortLines)

        for (i, firstLine) in enumerate(firstLines):
            secondLine = secondLines[i]
            firstCols = firstLine.split('\t')
            secondCols = secondLine.split('\t')
            for (j, firstCol) in enumerate(firstCols):
                secondCol = secondCols[j]
                try:
                    firstCol = float(firstCol)
                    secondCol = float(secondCol)
                    self.assertAlmostEqual(firstCol, secondCol, places=2)
                except ValueError:
                    self.assertEqual(firstCol, secondCol)

    def getFullModelTableFile(self):
        return os.path.join(self.pcssConfig["pcss_directory"], "data", "models", "human2008ModelTable.txt")

    def getLargeFastaFile(self):
        return os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")

    def getLargeDefinedFastaFile(self):
        return os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffDefinedInputFasta.txt")

    def getLinesToCompare(self, firstFile, secondFile, sortLines=False):
        firstReader = pcssTools.PcssFileReader(firstFile)
        secondReader = pcssTools.PcssFileReader(secondFile)
        if (sortLines):
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())
        else:
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())
        return [firstLines, secondLines]

class TestSequenceFeatures(PcssTest):

    def updateConfig(self, name, value):
        self.pcssConfig['name'] = value

    def readProteins(self):
        
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])

    def processException(self, peptide, errorAttribute, exceptionCode, function, *args):
        function(*args)

        self.assertEquals(peptide.getAttributeOutputString(errorAttribute), exceptionCode)
        self.assertRegexpMatches(peptide.getAttributeOutputString("peptide_errors"), exceptionCode)

    def getStringFeatureName(self):
        return "%s_string_feature" % self.name

    def getScoreFeatureName(self):
        return "%s_score_feature" % self.name

    def getPeptideNotFoundCode(self):
        return "peptide_error_%s_peptide_not_found" % self.name

    def getProteinMismatchCode(self):
        return "peptide_error_%s_protein_mismatch" % self.name

    def getBadCommandCode(self):
        return "peptide_error_%s_bad_command" % self.name

    def getBadCallCode(self):
        return "peptide_error_%s_bad_call" % self.name

    def getBadLineCode(self):
        return "peptide_error_%s_bad_line" % self.name

    def test_seq_feature_long(self):

        initialCwd = os.getcwd()
        self.setupLongRootDir()

        self.processNormalSeqFeature()
        time.sleep(5)
        shutil.rmtree(self.fileHandler.getSequenceFeatureDir(self.proteins[0].modbaseSequenceId))
        self.assertEquals(os.getcwd(), initialCwd)

    def test_seq_feature_short(self):
        self.processNormalSeqFeature()

    def processNormalSeqFeature(self):
        
        self.processResultFile()

        callString = self.getCallString()
        self.assertEquals(callString, self.seqData.getExpectedFullStringResult())
        
        peptide = self.proteins[0].peptides[17]
        self.assertEquals(peptide.attributes[self.getStringFeatureName()].getValueString(), self.seqData.stringFeatureValue)
        self.assertEquals(peptide.attributes[self.getScoreFeatureName()].getValueString(), self.seqData.scoreFeatureValue)

        self.proteins[0].setPeptide(pcssPeptide.PcssPeptide("FAKE", 1000, 1003, self.runner))
        self.processException(self.proteins[0].peptides[1000], self.getStringFeatureName(), self.getPeptideNotFoundCode(), self.processResultFile)
        
    def test_found_disopred_file(self):
        self.assertTrue(self.fileHandler.outputFileExists("76c3a409540532138c6b44bde9e4d248MDDRDENQ"))

    def test_residue_mismatch(self):
        
        self.pcssConfig["root_%s_dir" % self.name] = self.getErrorInputFile("residueMismatch")
        self.processException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getProteinMismatchCode(), self.processResultFile)

    def test_unexpected_call(self):

        self.pcssConfig["root_%s_dir" % self.name] = self.getErrorInputFile("badCall")
        self.processException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getBadCallCode(), self.processResultFile)

    def test_command_error(self):
        initialCwd = os.getcwd()
        self.setBadCommandData()
        self.assertRaises(pcssErrors.PcssGlobalException, self.processResultFile)
        self.assertEquals(os.getcwd(), initialCwd)

    def test_bad_line(self):
        self.pcssConfig["root_%s_dir" % self.name] = self.getErrorInputFile("badLine")
        self.processException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getBadLineCode(), self.processResultFile)


