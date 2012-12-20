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

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TestSvm(unittest.TestCase):
    
    def test_read_benchmark(self):
        br = pcssSvm.BenchmarkResults()
        br.readBenchmarkFile(self.pcssConfig['svm_benchmark_file'])

        st = br.getClosestScoreTuple(1)


    def getProtein(self, seqId, proteins):
        for protein in proteins:
            if (protein.modbaseSequenceId == seqId):
                return protein

    def getMaxFeatureNumber(self, proteinId, peptideId, reader):
        protein = self.getProtein(proteinId, reader.getProteins())
        peptide = protein.peptides[peptideId]
        svmLine = peptide.makeSvmFileLine()
        svmLineCols = svmLine.split(" ")
        lastFeatureNumber = svmLineCols[-1].split(":")[0]
        return lastFeatureNumber

    def test_write_svm_file(self):
        
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/annotationOutput.txt")
        modelFileName = self.pcssConfig["svm_model_file"]
        benchmarkFileName = self.pcssConfig["svm_benchmark_file"]
        appSvm = pcssSvm.ApplicationSvm(self.runner, modelFileName, benchmarkFileName)
        proteins = reader.getProteins()
        print "read %s proteins" % len(proteins)
        appSvm.setProteins(reader.getProteins())
        appSvm.writeClassificationFile()
        structureMaxFeatureNumber = self.getMaxFeatureNumber("ffb930a1b85cc26007aae5956ddf888dMEAFKKLR", 100, reader)
        self.assertEquals(structureMaxFeatureNumber, "192")
        
        nonStructureMaxFeatureNumber = self.getMaxFeatureNumber("ffc32c19c14fd6006e402bbf3a43c493MASTGYYA", 302, reader)
        self.assertEquals(nonStructureMaxFeatureNumber, "176")
        appSvm.classifySvm()
        
        appSvm.readResultFile()
        #except Exception as e:
        #print e
        #print e.msg

    def test_bad_svm_command(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/annotationOutput.txt")
        modelFileName = self.pcssConfig["svm_model_file"]
        benchmarkFileName = self.pcssConfig["svm_benchmark_file"]
        appSvm = pcssSvm.ApplicationSvm(self.runner, modelFileName, benchmarkFileName)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        appSvm.writeClassificationFile()
        self.runner.internalConfig["application_set_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], 
                                                                               "testInput/svmErrors/badCommandApplicationSet")
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.classifySvm)

        self.runner.internalConfig["application_set_file_name"] = "fake"
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.classifySvm)




    def test_bad_output_files(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/annotationOutput.txt")
        modelFileName = self.pcssConfig["svm_model_file"]
        benchmarkFileName = self.pcssConfig["svm_benchmark_file"]
        appSvm = pcssSvm.ApplicationSvm(self.runner, modelFileName, benchmarkFileName)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        appSvm.writeClassificationFile()
        
        self.runner.internalConfig["application_set_output_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], 
                                                                                      "testInput/svmErrors/countMismatchApplicationSet")
        try:
            appSvm.readResultFile()
        except Exception as e:
            print e.msg
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.readResultFile)
       
        self.runner.internalConfig["application_set_output_file_name"] = "fake"
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.readResultFile)

    def setUp(self):

        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"        

        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)


        self.runner = pcssTools.PcssRunner(self.pcssConfig)


    def test_invalid_benchmark_file(self):
        self.pcssConfig["svm_benchmark_file"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                             "testInput/svmErrors/invalidBenchmarkFile.txt")
        br = pcssSvm.BenchmarkResults()
        self.assertRaises(pcssErrors.PcssGlobalException, br.readBenchmarkFile, self.pcssConfig['svm_benchmark_file'])


    def test_nonstandard_aa(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmErrors/nonStandardAaAnnotation.txt")
        modelFileName = self.pcssConfig["svm_model_file"]
        benchmarkFileName = self.pcssConfig["svm_benchmark_file"]
        appSvm = pcssSvm.ApplicationSvm(self.runner, modelFileName, benchmarkFileName)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.writeClassificationFile)

    def test_invalid_svm_feature(self):
        self.runner.internalConfig["feature_order"][1] = "disorped_score_feature"
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/annotationOutput.txt")
        modelFileName = self.pcssConfig["svm_model_file"]
        benchmarkFileName = self.pcssConfig["svm_benchmark_file"]
        appSvm = pcssSvm.ApplicationSvm(self.runner, modelFileName, benchmarkFileName)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.writeClassificationFile)


if __name__ == '__main__':
    unittest.main()
