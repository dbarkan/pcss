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
import traceback
class PcssTest(unittest.TestCase):

    def setUp(self):
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"        
        #self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.pcssConfig = configobj.ConfigObj(configFile)
        oldRunDir = os.path.join(self.pcssConfig["home_test_directory"], "runs", self.pcssConfig["run_name"])

        if (os.path.exists(oldRunDir)):
            shutil.rmtree(oldRunDir)
        self.testExceptionOutputFh = open("testExceptionOutput.txt", "a")
        self.setupSpecificTest()

    def tearDown(self):
        self.testExceptionOutputFh.close()

    def getProtein(self, proteinId, proteins):
        for protein in proteins:
            if (protein.modbaseSequenceId == proteinId):
                return protein

    def clearErrorFiles(self):
        self.runner.pdh.runSubprocess(["rm", self.runner.pdh.getPcssErrorFile()], False)
        self.runner.pdh.runSubprocess(["rm", self.runner.pdh.getInternalErrorFile()], False)

    def getFullExpectedOutputFile(self, shortName):
        fullExpectedFileName = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "expectedOutput", 
                                            self.testName, "%s_expectedOutput.txt" % shortName)
        return fullExpectedFileName

    def createNewExpectedOutputFiles(self, observedFileName, shortExpectedFileName):
        expectedFileName = "%s_expectedOutput.txt" % shortExpectedFileName
        fullExpectedFileName = self.getFullExpectedOutputFile(shortExpectedFileName)

        outputDir = os.path.join(self.pcssConfig["home_test_directory"], "newExpectedOutputFiles")
        if (not os.path.exists(outputDir)):
            os.mkdir(outputDir)
            self.runner.sleepUntilDone(outputDir, predicate=pcssTools.pcssRunner.fileDoesNotExist)
        
        #copy observed file to temporary directory for holding new expected output files
        shutil.copy(observedFileName, os.path.join(outputDir, expectedFileName))

        #write command for copying file in temporary directory to expected output directory in test
        copyFh = open(os.path.join(outputDir, "copyFile.sh"), 'a')
        sourceFile = os.path.join(outputDir, expectedFileName)
        destinationFile = fullExpectedFileName
        copyFh.write("cp %s %s\n" % (sourceFile, destinationFile))

    def compareToExpectedOutput(self, observedFileName, shortExpectedFileName, sortLines=False, compareAlmostEqual=False):
        fullExpectedFileName = self.getFullExpectedOutputFile(shortExpectedFileName)
        #self.createNewExpectedOutputFiles(observedFileName, shortExpectedFileName)
        if (compareAlmostEqual):
            self.compareFilesAlmostEqual(observedFileName, fullExpectedFileName, sortLines)
        else:
            self.compareFiles(observedFileName, fullExpectedFileName, sortLines)

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
        return os.path.join(self.runner.internalConfig["pcss_directory"], "data", "models", "human2008ModelTable.txt")

    def getLargeFastaFile(self):
        return os.path.join(self.pcssConfig["user_pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")

    def getLargeDefinedFastaFile(self):
        return os.path.join(self.pcssConfig["user_pcss_directory"], "data", "inputSequences", "ffDefinedInputFasta.txt")

    def handleTestException(self, e):
        self.testExceptionOutputFh.write(e.exception.msg + "\n")
        tb = traceback.format_exc()
        self.testExceptionOutputFh.write(tb + "\n\n")

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

    def processFeatureException(self, peptide, errorAttribute, exceptionCode, function, *args):
        print "process exception args: %s" % str(*args)
        function(*args)

        self.assertEquals(peptide.getAttributeOutputString(errorAttribute), exceptionCode)
        self.assertRegexpMatches(peptide.getAttributeOutputString("peptide_errors"), exceptionCode)


class TestSequenceFeatures(PcssTest):

    def readProteins(self):
        
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.internalConfig = self.runner.internalConfig
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])

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

    def dtest_seq_feature_long(self):

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
        self.processFeatureException(self.proteins[0].peptides[1000], self.getStringFeatureName(), self.getPeptideNotFoundCode(), self.processResultFile)
        
    def test_found_disopred_file(self):
        self.assertTrue(self.fileHandler.outputFileExists("76c3a409540532138c6b44bde9e4d248MDDRDENQ"))

    def test_residue_mismatch(self):
        
        self.internalConfig["root_%s_dir" % self.name] = self.getErrorInputFile("residueMismatch")
        self.processFeatureException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getProteinMismatchCode(), self.processResultFile)

    def test_unexpected_call(self):

        self.internalConfig["root_%s_dir" % self.name] = self.getErrorInputFile("badCall")
        self.processFeatureException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getBadCallCode(), self.processResultFile)

    def test_command_error(self):
        initialCwd = os.getcwd()
        self.setBadCommandData()
        self.assertRaises(pcssErrors.PcssGlobalException, self.processResultFile)
        self.assertEquals(os.getcwd(), initialCwd)

    def test_bad_line(self):
        self.internalConfig["root_%s_dir" % self.name] = self.getErrorInputFile("badLine")
        self.processFeatureException(self.proteins[0].peptides.values()[0], self.getStringFeatureName(), self.getBadLineCode(), self.processResultFile)


