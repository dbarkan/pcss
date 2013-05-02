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

class TestSvm(pcssTests.PcssTest):
    
    def setupSpecificTest(self):
        self.errorDirName = "svmErrors"
        self.testName = "svm"

    def test_read_benchmark(self):
        br = pcssSvm.BenchmarkResults()
        br.readBenchmarkFile(self.pcssConfig['svm_benchmark_file'])
        st = br.getClosestScoreTuple(1)
          
        self.assertEquals(str(round(st.score, 3)), '0.998')

    def getMaxFeatureNumber(self, proteinId, peptideId, appSvm):
        protein = self.getProtein(proteinId, appSvm.getProteins())
        peptide = protein.peptides[peptideId]
        svmLine = peptide.makeSvmFileLine()
        svmLineCols = svmLine.split(" ")
        lastFeatureNumber = svmLineCols[-1].split(":")[0]
        return lastFeatureNumber

    def readTrainingAnnotationInputFile(self, fileName):
        self.setTrainingFileAttributes()
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        self.reader = pcssIO.AnnotationFileReader(self.runner)
        self.reader.readAnnotationFile(fileName)

    def readStandardTrainingAnnotationInputFile(self):
        self.readTrainingAnnotationInputFile("testInput/svmTrainingAnnotationInput.txt")

    def getSvmBenchmarker(self):
        self.readStandardTrainingAnnotationInputFile()
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        return benchmarker

    def readStandardSvmApplicationInputFile(self):
        self.readSvmApplicationInputFile("testInput/svmApplicationAnnotationInput.txt")

    def readSvmApplicationInputFile(self, fileName):
        self.setSvmApplicationFileAttributes()
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        self.reader = pcssIO.AnnotationFileReader(self.runner)
        self.reader.readAnnotationFile(fileName)

    def getApplicationSvm(self):
        self.readStandardSvmApplicationInputFile()
        appSvm = pcssSvm.ApplicationSvm(self.runner)
        proteins = self.reader.getProteins()
        appSvm.setProteins(self.reader.getProteins())
        return appSvm

    def getLeaveOneOutBenchmarker(self):

        self.readStandardTrainingAnnotationInputFile()
        benchmarker = pcssSvm.LeaveOneOutBenchmarker(self.runner)
        return benchmarker

    def test_write_svm_file(self):
        appSvm = self.getApplicationSvm()
        appSvm.writeClassificationFile()
        structureMaxFeatureNumber = self.getMaxFeatureNumber("ffb930a1b85cc26007aae5956ddf888dMEAFKKLR", 100, appSvm)
        self.assertEquals(structureMaxFeatureNumber, "192")

        nonStructureMaxFeatureNumber = self.getMaxFeatureNumber("ffc32c19c14fd6006e402bbf3a43c493MASTGYYA", 302, appSvm)
        self.assertEquals(nonStructureMaxFeatureNumber, "176")
        appSvm.classifySvm()

        appSvm.readResultFile()
        appSvm.addScoresToPeptides()
        protein = self.getProtein("ffb930a1b85cc26007aae5956ddf888dMEAFKKLR",  appSvm.getProteins())
        peptide = protein.peptides[100]
        self.assertEquals(peptide.getAttributeOutputString("svm_score"),  -1.6814754)

    def test_train_svm(self):
        benchmarker = self.getSvmBenchmarker()

        benchmarker.createTrainingAndTestSets(pcssTools.getAllPeptides(self.reader.getProteins(), False))
        self.assertEquals(len(benchmarker.benchmarkHandler.positiveTrainingSet), 61)
        self.assertEquals(len(benchmarker.benchmarkHandler.negativeTrainingSet), 61)
        self.assertEquals(len(benchmarker.benchmarkHandler.positiveTestSet), 6)
        self.assertEquals(len(benchmarker.benchmarkHandler.negativeTestSet), 73)
        
    def test_train_and_test_svm(self):
        benchmarker = self.getSvmBenchmarker()
        self.runner.internalConfig["make_random_test_set"] = False
        pstList = self.getSvmBenchmarkTestSetResult(benchmarker, pcssTools.getAllPeptides(self.reader.getProteins(), False))
        firstTuple = pstList.getBenchmarkTuple(0)
        self.assertEquals(float(firstTuple.score), -3.9027991)
        self.assertEquals(float(firstTuple.fpr), 0.0136986301369863)
        
    def test_leave_one_out(self):
        benchmarker = self.getLeaveOneOutBenchmarker()
        peptideSet = pcssTools.getAllPeptides(self.reader.getProteins(), False)[0:20]
        for i in range(len(peptideSet)):
            benchmarker.createTrainingAndTestSets(peptideSet)
            benchmarker.trainAndApplyModel() 
            benchmarker.readBenchmarkResults()
        benchmarker.processAllResults()
        self.compareToExpectedOutput(self.runner.pdh.getLeaveOneOutResultFileName(), "leaveOneOut")

    def test_leave_one_out_internal_count_error(self):

        benchmarker = self.getLeaveOneOutBenchmarker()
        peptideSet = pcssTools.getAllPeptides(self.reader.getProteins(), False)[0:20]
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            self.throwBadCountError(benchmarker, peptideSet)
        self.handleTestException(pge)

    def test_leave_one_out_test_set_error(self):

        benchmarker = self.getLeaveOneOutBenchmarker()
        peptideSet = pcssTools.getAllPeptides(self.reader.getProteins(), False)[0:20]
        benchmarker.createTrainingAndTestSets(peptideSet)
        fakeTestSvm = pcssSvm.TestSvm(self.runner)
        fakeTestSvm.setPeptides(pcssTools.getAllPeptides(self.reader.getProteins(), False)[21:23])
        benchmarker.testSvm = fakeTestSvm
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            benchmarker.trainAndApplyModel()
        self.handleTestException(pge)

    def throwBadCountError(self, benchmarker, peptideSet):
        badCount = len(peptideSet) + 1
        for i in range(badCount):
            benchmarker.createTrainingAndTestSets(peptideSet)
            benchmarker.trainAndApplyModel() 
            benchmarker.readBenchmarkResults()

    def getSvmBenchmarkTestSetResult(self, benchmarker, peptideList):
        benchmarker.createTrainingAndTestSets(peptideList)
        benchmarker.trainAndApplyModel() 
        benchmarker.testSvm.readResultFile()
        return benchmarker.testSvm.getTestSetResult()
        
    def test_multiple_iteration_normal_output(self):
                            
        benchmarker = self.getSvmBenchmarker()
        self.runner.internalConfig["make_random_test_set"] = False
        tsrt = pcssSvm.TestSetResultTracker(self.runner)

        peptideList = pcssTools.getAllPeptides(self.reader.getProteins(), False)
        testSetResult = self.getSvmBenchmarkTestSetResult(benchmarker, peptideList)
        tsrt.addTestSetResult(testSetResult)
        
        peptideList.reverse()
        testSetResult = self.getSvmBenchmarkTestSetResult(benchmarker, peptideList)
        tsrt.addTestSetResult(testSetResult)
        
        tsrt.finalize()
        tsrt.writeResultFile(self.runner.pdh.getFullBenchmarkResultFileName())
        self.compareToExpectedOutput(self.runner.pdh.getFullBenchmarkResultFileName(), "trainingSvm")
        
    def test_tracker_average(self):

        benchmarker = self.getSvmBenchmarker()
        tsrList = []
        for i in range(2):
            tsrList.append(self.getSvmBenchmarkTestSetResult(benchmarker, pcssTools.getAllPeptides(self.reader.getProteins(), False)))

        fprRates = []    
        tsrt = pcssSvm.TestSetResultTracker(self.runner)
        
        for tsr in tsrList:
            tprTuples = tsr.getIncrementedTprTuples()
            firstFprRate = tprTuples[0].fpr
            fprRates.append(firstFprRate)
            tsrt.addTestSetResult(tsr)

        average = float(sum(fprRates)) / float(len(fprRates))
        tsrt.finalize()
        firstTuple = tsrt.getBenchmarkTuple(0)
        self.assertAlmostEqual(average, firstTuple.fpr)

    def test_test_sets_count_positive_mismatch(self):
        runner = pcssTools.PcssRunner(self.pcssConfig)
        tsrt = self.countMismatch(runner.getPositiveKeyword(), 8)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            tsrt.finalize()
        self.handleTestException(e)
        self.assertTrue(runner.getPositiveKeyword().lower() in e.exception.msg)

    def test_test_sets_count_negative_mismatch(self):
        runner = pcssTools.PcssRunner(self.pcssConfig)
        tsrt = self.countMismatch(runner.getNegativeKeyword(), 1)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            tsrt.finalize()
        self.handleTestException(e)
        self.assertTrue(runner.getNegativeKeyword().lower() in e.exception.msg)

    def countMismatch(self, keyword, peptidesToRemoveCount):
        
        benchmarker = self.getSvmBenchmarker()
        self.runner.internalConfig["make_random_test_set"] = False
        peptideList = pcssTools.getAllPeptides(self.reader.getProteins(), False)
        tsrt = pcssSvm.TestSetResultTracker(self.runner)

        firstTestSetResult = self.getSvmBenchmarkTestSetResult(benchmarker, peptideList)

        positivePeptidesRemoved = 0
        for protein in self.reader.getProteins():
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

        peptideList = pcssTools.getAllPeptides(self.reader.getProteins(), False)
        secondTestSetResult = self.getSvmBenchmarkTestSetResult(benchmarker, peptideList)
        tsrt.addTestSetResult(firstTestSetResult)
        tsrt.addTestSetResult(secondTestSetResult)
        return tsrt
            
    def test_bad_status(self):
        self.readStandardTrainingAnnotationInputFile()
        
        peptide = self.reader.getProteins()[0].peptides.values()[0]
        peptide.addStringAttribute("status", "fake")
        svm = pcssSvm.TrainingSvm(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            svm.setPeptides(pcssTools.getAllPeptides(self.reader.getProteins(), False))
        self.handleTestException(pge)

    def test_bad_svm_training_command(self):
        self.readStandardTrainingAnnotationInputFile()
        svm = pcssSvm.TrainingSvm(self.runner)
        self.runner.internalConfig["training_set_file_name"] = self.getErrorInputFile("badCommandTrainingSet")

        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            svm.trainModel()
        self.handleTestException(pge)

        self.runner.internalConfig["training_set_file_name"] = "fake"
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            svm.trainModel()
        self.handleTestException(pge)

    def test_more_positives_than_negatives(self):
        self.readTrainingAnnotationInputFile(self.getErrorInputFile("trainingMorePositives.txt"))
        peptides = pcssTools.getAllPeptides(self.reader.getProteins(), False)
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            benchmarker.createTrainingAndTestSets(peptides)
        self.handleTestException(pge)

    def test_no_test_set_positives(self):
        try:
            self.readTrainingAnnotationInputFile(self.getErrorInputFile("trainingNoTestSetPositives.txt"))
        except pcssErrors.PcssGlobalException as e:
            print e.msg
        peptides = pcssTools.getAllPeptides(self.reader.getProteins(), False)
        benchmarker = pcssSvm.SvmBenchmarker(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            benchmarker.createTrainingAndTestSets(peptides)
        self.handleTestException(pge)

    def test_bad_svm_app_command(self):

        appSvm = self.getApplicationSvm()
        appSvm.writeClassificationFile()
        self.runner.internalConfig["application_set_file_name"] = self.getErrorInputFile("badCommandApplicationSet")

        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.classifySvm()
        self.handleTestException(pge)

        self.runner.internalConfig["application_set_file_name"] = "fake"

        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.classifySvm()
        self.handleTestException(pge)

    def test_missing_training_model(self):
        self.readStandardTrainingAnnotationInputFile()

        testSvm = pcssSvm.TestSvm(self.runner)
        testSvm.setProteins(self.reader.getProteins())
        testSvm.writeClassificationFile()

        self.runner.internalConfig["training_new_model_name"] = "fake"
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            testSvm.classifySvm()
        self.handleTestException(pge)

    def test_bad_app_output_files(self):
        appSvm = self.getApplicationSvm()
        appSvm.writeClassificationFile()
        
        self.runner.internalConfig["application_set_output_file_name"] = self.getErrorInputFile("countMismatchApplicationSet")
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.readResultFile()
        self.handleTestException(pge)
        
        self.runner.internalConfig["application_set_output_file_name"] = "fake"
        
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.readResultFile()
        self.handleTestException(pge)

    
    def processInvalidBenchmarkFile(self, fileName):
        br = pcssSvm.BenchmarkResults()        
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            br.readBenchmarkFile(fileName)
        self.handleTestException(pge)

    def test_all_invalid_benchmark_files(self):
        self.processInvalidBenchmarkFile(self.getErrorInputFile("invalidBenchmarkFile.txt"))
        self.processInvalidBenchmarkFile(self.getErrorInputFile("benchmarkFileNoFirstLine.txt"))
        self.processInvalidBenchmarkFile(self.getErrorInputFile("benchmarkFileNoLastLine.txt"))
        
    def test_nonstandard_aa(self):
        self.readSvmApplicationInputFile(self.getErrorInputFile("nonStandardAaAnnotation.txt"))
        appSvm = pcssSvm.ApplicationSvm(self.runner)
        proteins = self.reader.getProteins()
        appSvm.setProteins(self.reader.getProteins())
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.writeClassificationFile()
        self.handleTestException(pge)

    def test_invalid_svm_feature(self):

        appSvm = self.getApplicationSvm()
        self.runner.internalConfig["feature_order"][1] = "disorped_score_feature"
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            appSvm.writeClassificationFile()
        self.handleTestException(pge)

if __name__ == '__main__':
    unittest.main()
