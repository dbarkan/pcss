import unittest
import time
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
#import objgraph

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestRunner(unittest.TestCase):

    def setUp(self):
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)

    def clearErrorFiles(self):
        self.runner.pdh.runSubprocess(["rm", self.runner.pdh.getPcssErrorFile()], False)
        self.runner.pdh.runSubprocess(["rm", self.runner.pdh.getInternalErrorFile()], False)

    def test_annotation_pcss_error(self):
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.throwPcssException()

    def test_svm_application_feature_pcss_error(self):
        self.runner = pcssTools.SvmApplicationFeatureRunner(self.pcssConfig)
        self.throwPcssException()
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

    def test_svm_application_input_pcss_error(self):
        self.pcssConfig["input_annotation_file_name"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                                     "testInput/ioErrors/missingColumnsFile.txt")
        self.runner = pcssTools.SvmApplicationInputRunner(self.pcssConfig)
        
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

    def test_svm_training_pcss_error(self):
        self.pcssConfig["input_annotation_file_name"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                                     "testInput/ioErrors/missingColumnsFile.txt")
        self.runner = pcssTools.TrainingBenchmarkRunner(self.pcssConfig)
        
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)
        
    
    def throwPcssException(self):

        self.pcssConfig["peptide_importer_type"] = "defined"
        self.pcssConfig['fasta_file'] = "testInput/ioErrors/peptideMismatchFasta.txt"
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        

    def test_annotation_internal_error(self):
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.throwInternalException()
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

    def test_svm_application_feature_internal_error(self):
        self.runner = pcssTools.SvmApplicationFeatureRunner(self.pcssConfig)
        self.throwInternalException()
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

    def test_svm_application_input_internal_error(self):
        
        self.runner = pcssTools.SvmApplicationInputRunner(self.pcssConfig)
        self.pcssConfig['fasta_file'] = "fake"
        print "RUN SVM INPUT"

        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

    def test_svm_training_input_internal_error(self):

        self.runner = pcssTools.TrainingBenchmarkRunner(self.pcssConfig)
        self.pcssConfig["jackknife_fraction"] = "fake"
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.assertRaises(pcssErrors.ErrorExistsException, self.runner.execute)

        
    def throwInternalException(self):
        self.pcssConfig["peptide_importer_type"] = "defined"
        self.pcssConfig['fasta_file'] = "fake"
        self.clearErrorFiles()
        self.runner.execute()
        self.assertTrue(os.path.exists(self.runner.pdh.getInternalErrorFile()))

    def test_config_error(self):
        configFile = "testInput/ioErrors/pcssConfigMissingFile.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        print "testing config error"
        obj = configobj.ConfigObj(configFile, configspec=configSpecFile)
        self.assertRaises(pcssErrors.PcssGlobalException, pcssTools.SvmApplicationFeatureRunner, obj)

    def test_svm_application_features(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "svmApplicationFileAttributes.txt")
        self.setLargeFastaFile()
        self.runner = pcssTools.SvmApplicationFeatureRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareFilesAlmostEqual(self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]),
                          self.getExpectedOutputFile("svmApplication"), True)

    def test_svm_application_input(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "svmApplicationFileAttributes.txt")
        self.setLargeFastaFile()
        self.runner = pcssTools.SvmApplicationInputRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareFilesAlmostEqual(self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]),
                                     self.getExpectedOutputFile("svmApplication"), True)

    def setLargeFastaFile(self):
        #self.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFastaLong.txt")
        self.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "inputSequences", "ffSequencesFasta.txt")
        self.pcssConfig['model_table_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "models", "human2008ModelTable.txt")

    def setLargeDefinedFastaFile(self):
        self.pcssConfig['model_table_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "models", "human2008ModelTable.txt")
        self.pcssConfig["fasta_file"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "ffDefinedInputFasta.txt")

    def test_training_annotation_runner(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.pcssConfig["peptide_importer_type"] = "defined"
        self.setLargeDefinedFastaFile()
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareFilesAlmostEqual(self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]),
                                     self.getExpectedOutputFile("trainingAnnotation"))

    def test_training_svm_runner(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.pcssConfig["input_annotation_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "svmTrainingAnnotationInput.txt")
        self.runner = pcssTools.TrainingBenchmarkRunner(self.pcssConfig)
        #self.runner.internalConfig["make_random_test_set"] = False
        self.clearErrorFiles()
        self.runner.internalConfig["make_random_test_set"] = False
        self.runner.execute()

    def test_create_model(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.pcssConfig["input_annotation_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "svmTrainingAnnotationInput.txt")
        self.runner = pcssTools.CompleteSvmRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.compareFiles(self.runner.pdh.getUserModelFileName(), self.getExpectedOutputFile("userModel"))

    def dtest_leave_one_out_runner(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.pcssConfig["input_annotation_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "svmTrainingAnnotationInput.txt")
        self.runner = pcssTools.LeaveOneOutBenchmarkRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.compareFiles(self.runner.pdh.getLeaveOneOutResultFileName(), self.getExpectedOutputFile("trainingLoo"))
    def test_annotation_runner(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "annotationFileAttributes.txt")
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.setLargeFastaFile()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareFilesAlmostEqual(self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]),
                                     self.getExpectedOutputFile("annotation"))

    def test_create_feature_error(self):
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        self.pcssConfig['model_table_file'] = os.path.join(self.pcssConfig["pcss_directory"], "data", "models", "fakeIdModelTable.txt")
        self.runner.execute()

    def compareFilesAlmostEqual(self, firstFile, secondFile, sortLines=False):
        firstReader = pcssTools.PcssFileReader(firstFile)
        secondReader = pcssTools.PcssFileReader(secondFile)
        if (sortLines):
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())
        else:
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())

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
                    #print "first col %s and second col %s are equal" % (firstCol, secondCol)
                except ValueError:
                    self.assertEqual(firstCol, secondCol)

    def compareFiles(self, firstFile, secondFile, sortLines=False):
        firstReader = pcssTools.PcssFileReader(firstFile)
        secondReader = pcssTools.PcssFileReader(secondFile)
        if (sortLines):
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())
        else:
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())

        for (i, firstLine) in enumerate(firstLines):

            secondLine = secondLines[i]
            self.assertEquals(firstLine, secondLine)



    def getExpectedOutputFile(self, filePrefix):
        return os.path.join(self.pcssConfig["home_test_directory"], "testInput", "expectedOutput", "%s_expectedOutput.txt" % filePrefix)


if __name__ == '__main__':
    unittest.main()
