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
        try:
            reader = pcssIO.AnnotationFileReader(self.runner)
            reader.readAnnotationFile("testInput/svmApplicationAnnotationInput.txt")

            appSvm = pcssSvm.ApplicationSvm(self.runner)
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
            appSvm.addScoresToPeptides()
            protein = self.getProtein("ffb930a1b85cc26007aae5956ddf888dMEAFKKLR",  reader.getProteins())
            peptide = protein.peptides[100]
            self.assertEquals(peptide.getAttributeOutputString("svm_score"),  -1.6814754)
        except pcssErrors.PcssGlobalException as e:
            print e.msg
    def test_train_svm(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
        self.assertEquals(len(benchmarker.benchmarkHandler.positiveTrainingSet), 61)
        self.assertEquals(len(benchmarker.benchmarkHandler.negativeTrainingSet), 61)
        self.assertEquals(len(benchmarker.benchmarkHandler.positiveTestSet), 6)
        self.assertEquals(len(benchmarker.benchmarkHandler.negativeTestSet), 73)
        
    def test_train_and_test_svm(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        self.runner.internalConfig["make_random_test_set"] = False
        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
        benchmarker.trainAndApplyModel() 
        benchmarker.readBenchmarkResults()
        pstList = benchmarker.testSvm.getTestSetResult()
        firstTuple = pstList.getBenchmarkTuple(0)
        self.assertEquals(float(firstTuple.score), -4.8898455)
        self.assertEquals(float(firstTuple.fpr), 0.0136986301369863)
        
    def test_leave_one_out(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.LeaveOneOutBenchmarker(self.runner)
        peptideSet = pcssTools.getAllPeptides(reader.getProteins(), False)[0:20]
        for i in range(len(peptideSet)):
            benchmarker.createTrainingAndTestSets(peptideSet)
            benchmarker.trainAndApplyModel() 
            benchmarker.readBenchmarkResults()
        benchmarker.processAllResults()
        expectedFile = os.path.join(self.runner.pcssConfig["home_test_directory"], "testInput", "expectedOutput", "leaveOneOut_expectedOutput.txt")
        observedOutputFile = self.runner.pdh.getLeaveOneOutResultFileName()
        self.compareFiles(expectedFile, observedOutputFile)

    def test_leave_one_out_internal_count_error(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.LeaveOneOutBenchmarker(self.runner)
        peptideSet = pcssTools.getAllPeptides(reader.getProteins(), False)[0:20]
        
        self.assertRaises(pcssErrors.PcssGlobalException, self.throwBadCountError, benchmarker, peptideSet)

    def test_leave_one_out_test_set_error(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.LeaveOneOutBenchmarker(self.runner)
        peptideSet = pcssTools.getAllPeptides(reader.getProteins(), False)[0:20]
        benchmarker.createTrainingAndTestSets(peptideSet)
        fakeTestSvm = pcssSvm.TestSvm(self.runner)
        fakeTestSvm.setPeptides(pcssTools.getAllPeptides(reader.getProteins(), False)[21:23])
        benchmarker.testSvm = fakeTestSvm
        self.assertRaises(pcssErrors.PcssGlobalException, benchmarker.trainAndApplyModel)

    def throwBadCountError(self, benchmarker, peptideSet):
        badCount = len(peptideSet) + 1
        for i in range(badCount):
            benchmarker.createTrainingAndTestSets(peptideSet)
            benchmarker.trainAndApplyModel() 
            benchmarker.readBenchmarkResults()


    def test_multiple_iteration_normal_output(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        tsrt = pcssSvm.TestSetResultTracker(self.runner)
        reader = pcssIO.AnnotationFileReader(self.runner)
        self.runner.internalConfig["make_random_test_set"] = False
        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        
        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
        benchmarker.trainAndApplyModel() 
        benchmarker.testSvm.readResultFile()
        

        tsrt.addTestSetResult(benchmarker.testSvm.getTestSetResult())
        
        newPeptideList = pcssTools.getAllPeptides(reader.getProteins(), False)
        newPeptideList.reverse()
        
        benchmarker.createTrainingAndTestSets(newPeptideList)
        benchmarker.trainAndApplyModel() 
        benchmarker.testSvm.readResultFile()
        
        tsrt.addTestSetResult(benchmarker.testSvm.getTestSetResult())
        
        tsrt.finalize()
        tsrt.writeResultFile(self.runner.pdh.getFullBenchmarkResultFileName())
        expectedOutputFile = os.path.join(self.runner.pcssConfig["home_test_directory"], "testInput", "expectedOutput", "trainingSvm_expectedOutput.txt")
        observedOutputFile = self.runner.pdh.getFullBenchmarkResultFileName()
        self.compareFiles(expectedOutputFile, observedOutputFile)
        

    def test_tracker_average(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        
        tsrList = []
        for i in range(2):
            benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
            benchmarker.trainAndApplyModel() 
            benchmarker.testSvm.readResultFile()
            tsrList.append(benchmarker.testSvm.getTestSetResult())
        fprRates = []    
        tsrt = pcssSvm.TestSetResultTracker(self.runner)
        
        for tsr in tsrList:
            tprTuples = tsr.getIncrementedTprTuples()
            firstFprRate = tprTuples[0].fpr
            fprRates.append(firstFprRate)
            tsrt.addTestSetResult(tsr)
        print "LIST OF FPR RATES: %s" % fprRates
        average = float(sum(fprRates)) / float(len(fprRates))
        tsrt.finalize()
        firstTuple = tsrt.getBenchmarkTuple(0)
        self.assertAlmostEqual(average, firstTuple.fpr)

    def test_test_sets_count_positive_mismatch(self):
        tsrt = self.countMismatch(self.runner.getPositiveKeyword(), 8)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            tsrt.finalize()
        print "GOT EXCPETION %s" % e.exception.msg
        self.assertTrue(self.runner.getPositiveKeyword().lower() in e.exception.msg)

    def test_test_sets_count_negative_mismatch(self):
        tsrt = self.countMismatch(self.runner.getNegativeKeyword(), 1)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            tsrt.finalize()
        print "GOT EXCPETION %s" % e.exception.msg
        self.assertTrue(self.runner.getNegativeKeyword().lower() in e.exception.msg)

    def countMismatch(self, keyword, peptidesToRemoveCount):
        
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
                            
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        self.runner.internalConfig["make_random_test_set"] = False
        
        tsrList = []
        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
        benchmarker.trainAndApplyModel() 
        benchmarker.testSvm.readResultFile()
        firstTestSetResult = benchmarker.testSvm.getTestSetResult()
        positivePeptidesRemoved = 0
        for protein in reader.getProteins():
            peptides = protein.peptides

            peptideIter = peptides.iteritems()
            for peptidePosition, peptide in list(peptideIter):
                if (peptide.getAttributeOutputString("status") == keyword):
                    peptides.pop(peptidePosition)
                    positivePeptidesRemoved += 1
                    if (positivePeptidesRemoved > peptidesToRemoveCount):
                        break
            if (positivePeptidesRemoved > peptidesToRemoveCount):
                break

        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(reader.getProteins(), False))
        benchmarker.trainAndApplyModel() 
        benchmarker.testSvm.readResultFile()
        secondTestSetResult = benchmarker.testSvm.getTestSetResult()
                
        tsrt = pcssSvm.TestSetResultTracker(self.runner)
        tsrt.addTestSetResult(firstTestSetResult)
        tsrt.addTestSetResult(secondTestSetResult)
        return tsrt
            
    def test_bad_status(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
        peptide = reader.getProteins()[0].peptides.values()[0]
        peptide.addStringAttribute("status", "fake")
        svm = pcssSvm.TrainingSvm(self.runner)
        self.assertRaises(pcssErrors.PcssGlobalException, svm.setPeptides, pcssTools.getAllPeptides(reader.getProteins(), False))

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


    def test_bad_svm_training_command(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
        svm = pcssSvm.TrainingSvm(self.runner)
        self.runner.internalConfig["training_set_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], 
                                                                            "testInput/svmErrors/badCommandTrainingSet")
        self.assertRaises(pcssErrors.PcssGlobalException, svm.trainModel)
        self.runner.internalConfig["training_set_file_name"] = "fake"
            
        self.assertRaises(pcssErrors.PcssGlobalException, svm.trainModel)

    def test_more_positives_than_negatives(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmErrors/trainingMorePositives.txt")
        peptides = pcssTools.getAllPeptides(reader.getProteins(), False)
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)

        self.assertRaises(pcssErrors.PcssGlobalException, benchmarker.createTrainingAndTestSets, peptides)

    def test_no_test_set_positives(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmErrors/trainingNoTestSetPositives.txt")
        peptides = pcssTools.getAllPeptides(reader.getProteins(), False)
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        self.assertRaises(pcssErrors.PcssGlobalException, benchmarker.createTrainingAndTestSets, peptides)

    def test_bad_svm_app_command(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmApplicationAnnotationInput.txt")
        appSvm = pcssSvm.ApplicationSvm(self.runner)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        appSvm.writeClassificationFile()
        self.runner.internalConfig["application_set_file_name"] = os.path.join(self.pcssConfig["home_test_directory"], 
                                                                               "testInput/svmErrors/badCommandApplicationSet")
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.classifySvm)

        self.runner.internalConfig["application_set_file_name"] = "fake"
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.classifySvm)

    def test_missing_training_model(self):
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "trainingFileAttributes.txt")
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile("testInput/svmTrainingAnnotationInput.txt")
        testSvm = pcssSvm.TestSvm(self.runner)
        proteins = reader.getProteins()
        testSvm.setProteins(reader.getProteins())
        testSvm.writeClassificationFile()

        self.runner.internalConfig["training_new_model_name"] = "fake"
        self.assertRaises(pcssErrors.PcssGlobalException, testSvm.classifySvm)

    def test_bad_app_output_files(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmApplicationAnnotationInput.txt")
        appSvm = pcssSvm.ApplicationSvm(self.runner)
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
        self.pcssConfig["attribute_file_name"] = os.path.join(self.pcssConfig["pcss_directory"], "data", "context", "svmApplicationFileAttributes.txt")

        self.runner = pcssTools.PcssRunner(self.pcssConfig)


    def test_invalid_benchmark_file(self):
        self.pcssConfig["svm_benchmark_file"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                             "testInput/svmErrors/invalidBenchmarkFile.txt")
        br = pcssSvm.BenchmarkResults()
        self.assertRaises(pcssErrors.PcssGlobalException, br.readBenchmarkFile, self.pcssConfig['svm_benchmark_file'])
        
        self.pcssConfig["svm_benchmark_file"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                             "testInput/svmErrors/benchmarkFileNoFirstLine.txt")
        br = pcssSvm.BenchmarkResults()

        self.assertRaises(pcssErrors.PcssGlobalException, br.readBenchmarkFile, self.pcssConfig['svm_benchmark_file'])
        
        self.pcssConfig["svm_benchmark_file"] = os.path.join(self.pcssConfig["home_test_directory"],
                                                             "testInput/svmErrors/benchmarkFileNoLastLine.txt")
        br = pcssSvm.BenchmarkResults()
        self.assertRaises(pcssErrors.PcssGlobalException, br.readBenchmarkFile, self.pcssConfig['svm_benchmark_file'])


    def test_nonstandard_aa(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmErrors/nonStandardAaAnnotation.txt")
        appSvm = pcssSvm.ApplicationSvm(self.runner)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.writeClassificationFile)

    def test_invalid_svm_feature(self):
        self.runner.internalConfig["feature_order"][1] = "disorped_score_feature"
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmApplicationAnnotationInput.txt")
        appSvm = pcssSvm.ApplicationSvm(self.runner)
        proteins = reader.getProteins()
        appSvm.setProteins(reader.getProteins())
        self.assertRaises(pcssErrors.PcssGlobalException, appSvm.writeClassificationFile)


if __name__ == '__main__':
    unittest.main()
