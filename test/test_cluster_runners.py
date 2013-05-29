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
import shutil
import sys
import pcssTests
import pcssSvm
import pcssIO
from validate import Validator

class TestClusterRunner(pcssTests.PcssTest):
    def setupSpecificTest(self):
        self.testName = "clusterRunner"

    def getExpectedOutputFile(self, fileName):
        return os.path.join(self.pcssConfig["pcss_directory"], "test", "testInput", "expectedOutput", fileName)

    def dtest_prepare_svm_application_server_runner(self):
        userConfig = configobj.ConfigObj(os.path.join("testConfig", "testSvmApplicationServerConfig.txt"))

        jobDirectory = os.path.join("/trombone1/home/dbarkan/pcss/", "test", "runs", "develop")
        if (not os.path.exists(jobDirectory)):
            os.mkdir(jobDirectory)

        userConfig["job_directory"] = jobDirectory #also set by job class -- self.directory
        userConfig["run_name"] = "develop"         #equivalent to self.name
        print "job directory: %s" % jobDirectory
        self.runner = pcssTools.PrepareSvmApplicationServerRunner(userConfig)

        self.copyAnnotationServerInput(jobDirectory)
        self.runner.execute()
        shellScript =  self.runner.getClusterShellScript()
        outputFh = open(self.runner.pdh.getFullOutputFile("testShellScript"), "w")
        outputFh.write(shellScript)
        outputFh.close()
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("seqBatchList/0/seqBatchNodeClusterConfig.txt"), "svmApplicationServerSeqBatchConfig")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("testShellScript"), "svmApplicationServerSubmit.sh")

    def copyAnnotationServerInput(self, runDirectory):
        sourceDirectory = os.path.join(self.pcssConfig["home_test_directory"], "testInput")
        
        shutil.copy(os.path.join(sourceDirectory, "ffSequencesFasta.txt"), os.path.join(runDirectory, self.runner.internalConfig["server_input_fasta_file_name"]))
        shutil.copy(os.path.join(sourceDirectory, "peptideRulesFile"), os.path.join(runDirectory, self.runner.internalConfig["server_input_rules_file_name"]))

    def test_prepare_svm_application_cluster_runner(self):
        self.pcssConfig["fasta_file"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "ffSequencesFasta.txt")
        self.runner = pcssTools.PrepareSvmApplicationClusterRunner(self.pcssConfig)
        self.clearErrorFiles()

        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("svmApplicationSubmit.sh"), "svmApplicationSubmit.sh")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("seqBatchList/0/inputFastaFile.txt"), "svmApplicationSeqBatchFasta")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("seqBatchList/0/seqBatchNodeClusterConfig.txt"), "svmApplicationSeqBatchConfig")
    
    def dtest_finalize_svm_application_cluster_runner(self):
        self.pcssConfig["fasta_file"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "ffSequencesFasta.txt")
        self.runner = pcssTools.FinalizeApplicationClusterRunner(self.pcssConfig)
        self.copySeqBatchToRunDir()
        self.clearErrorFiles()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("annotationOutput.txt"), "finalizeSvmApplication", sortLines=True)
        
    def copySeqBatchToRunDir(self):

        destinationDir = self.runner.pdh.getFullOutputFile("seqBatchList")
        if (os.path.exists(destinationDir)):
            shutil.rmtree(destinationDir)
        sourceDir = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "cluster", "seqBatchList")
        shutil.copytree(sourceDir, destinationDir)

    def dtest_prepare_training_annotation_cluster_runner(self):
        
        self.pcssConfig["peptide_importer_type"] = "defined"
        self.pcssConfig["fasta_file"] = os.path.join(self.pcssConfig["home_test_directory"], "testInput", "ffDefinedInputFasta.txt")
        self.runner = pcssTools.PrepareTrainingAnnotationClusterRunner(self.pcssConfig)
        self.clearErrorFiles()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("trainingAnnotationSubmit.sh"), "trainingAnnotationSubmit.sh")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("seqBatchList/0/inputFastaFile.txt"), "trainingAnnotationSeqBatchFasta")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("seqBatchList/0/seqBatchNodeClusterConfig.txt"), "trainingAnnotationSeqBatchConfig")
        
    def dtest_prepare_training_benchmark_cluster_runner(self):

        self.runner = pcssTools.PrepareTrainingBenchmarkClusterRunner(self.pcssConfig)
        self.clearErrorFiles()
        
        tempFh = open(self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]), 'w')
        tempFh.close()
        self.runner.execute()
        self.assertFalse(os.path.exists(self.runner.pdh.getPcssErrorFile()))
        self.assertFalse(os.path.exists(self.runner.pdh.getInternalErrorFile()))
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("trainingBenchmarkSubmit.sh"), "trainingBenchmarkSubmit.sh")
        self.compareToExpectedOutput(self.runner.pdh.getFullOutputFile("trainingBenchmarkNodeClusterConfig.txt"), "trainingBenchmarkNodeConfig")


if __name__ == '__main__':
    unittest.main()
