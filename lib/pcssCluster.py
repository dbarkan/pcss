import sys
import os
from Bio import SeqIO
import itertools
import pcssIO
import pcssModels
import pcssTools
import copy
import logging
import pcssPeptide
import pcssErrors
import pcssFeatures
import pcssFeatureHandlers
log = logging.getLogger("pcssPeptide")

class SeqDivider:
    def __init__(self, pcssRunner):

        self.pdh = pcssRunner.pdh
        self.pcssRunner = pcssRunner
        self.subDirectories = []
        self.csg = ClusterScriptGenerator(pcssRunner)

    def getFastaGroupList(self, fastaFile):
        fh = open(fastaFile, 'r')
        fastaList = list(SeqIO.FastaIO.FastaIterator(fh))
        seqBatchSize = int(self.pcssRunner.internalConfig["seq_batch_size"])

        print "seq batch size: %s" % seqBatchSize
        seqGroupList = [fastaList[i:i+seqBatchSize] for i in range(0, len(fastaList), seqBatchSize)]
        fh.close()
        self.seqBatchCount = len(seqGroupList)
        return seqGroupList

    def getSequenceTaskListString(self):
        taskList = []
        for i in range(self.seqBatchCount):
            taskList.append(str(i))
        return " ".join(taskList)
 
    def mergeSvmApplicationResults(self):
        
        fastaFile = self.pcssRunner.pcssConfig['fasta_file']
        seqGroupList = self.getFastaGroupList(fastaFile)
        print "merging application results"
        seqBatchDirectory = self.pcssRunner.getSeqBatchDirectory()
        allProteins = []
        for (i, nextGroup) in enumerate(seqGroupList):
            subDirName = self.getSeqBatchSubDirectoryName(i)
            subOutputFile = os.path.join(subDirName, self.pcssRunner.internalConfig["annotation_output_file"])
            reader = pcssIO.AnnotationFileReader(self.pcssRunner)
            reader.readAnnotationFile(subOutputFile)
            proteins = reader.getProteins()
            allProteins += proteins
        
        afw = pcssIO.AnnotationFileWriter(self.pcssRunner)
        afw.writeAllOutput(allProteins)

    def mergeTrainingAnnotationResults(self):
        
        fastaFile = self.pcssRunner.pcssConfig['fasta_file']
        seqGroupList = self.getFastaGroupList(fastaFile)

        seqBatchDirectory = self.pcssRunner.getSeqBatchDirectory()
        allProteins = []
        for (i, nextGroup) in enumerate(seqGroupList):
            subDirName = self.getSeqBatchSubDirectoryName(i)
            subOutputFile = os.path.join(subDirName, self.pcssRunner.internalConfig["annotation_output_file"])
            reader = pcssIO.AnnotationFileReader(self.pcssRunner)
            reader.readAnnotationFile(subOutputFile)
            proteins = reader.getProteins()
            allProteins += proteins
        
        afw = pcssIO.AnnotationFileWriter(self.pcssRunner)
        afw.writeAllOutput(allProteins)


    #when we get here, we have fasta file and unannotated sequences. 
    def divideSeqsFromFasta(self):
        fastaFile = self.pcssRunner.pcssConfig['fasta_file']
        seqGroupList = self.getFastaGroupList(fastaFile)

        seqBatchDirectory = self.pcssRunner.getSeqBatchDirectory()
        
        for (i, nextGroup) in enumerate(seqGroupList):
            self.createSequenceSubDirectory(i)

            self.makeSeqBatchFile(i, nextGroup)
        
            self.makeParameterFile(i)

    def makeParameterFile(self, i):
        pcssCopy = copy.deepcopy(self.pcssRunner.pcssConfig)

        if "run_name" in pcssCopy:
            del pcssCopy["run_name"]
        if "run_directory" in pcssCopy:
            del pcssCopy["run_directory"]

        pcssCopy["fasta_file"] =  self.getSeqBatchFileName(i)
        pcssCopy["run_name"] = str(i)
        pcssCopy["run_directory"] = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["seq_batch_directory"])
        pcssCopy["using_web_server"] = False

        subDirectoryName = self.getSeqBatchSubDirectoryName(i)
        pcssCopy.filename = os.path.join(subDirectoryName, self.pcssRunner.internalConfig["seq_batch_parameter_file_name"])
        pcssCopy.write()
        
    def getSeqBatchSubDirectoryName(self, index):
        seqBatchDirectory = self.pcssRunner.internalConfig["seq_batch_directory"]
        return self.pcssRunner.pdh.getFullOutputFile(os.path.join(seqBatchDirectory, str(index)))
        
    def getSeqBatchDirectory(self):
        return self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["seq_batch_directory"])

    def createSequenceSubDirectory(self, i):
        subDirectoryName = self.getSeqBatchSubDirectoryName(i)
        if (not os.path.exists(subDirectoryName)):
            os.mkdir(subDirectoryName)

    def getSeqBatchFileName(self, i):
        subDirectoryName = self.getSeqBatchSubDirectoryName(i)
        fileName = os.path.join(subDirectoryName, "%s" % self.pcssRunner.internalConfig["seq_batch_input_fasta_file_name"])
        return fileName
        
    def makeSeqBatchFile(self, i, seqBatch):
        fileName = self.getSeqBatchFileName(i)
        fh = open(fileName, 'w')
        for seqRecord in seqBatch:
            fh.write(">%s\n" % seqRecord.id)
            fh.write("%s\n" % str(seqRecord.seq))
        fh.close()

    def makeFullSvmApplicationSgeScript(self):
        nodeScriptName = self.pcssRunner.internalConfig["svm_application_node_script"]
        shellScriptName = self.pcssRunner.internalConfig["svm_application_shell_script"]
        outputFileName = self.pcssRunner.internalConfig["svm_application_stdout_file"]
        self.makeFullAnnotationSgeScript(nodeScriptName, shellScriptName, outputFileName)

    def makeFullTrainingAnnotationSgeScript(self):
        nodeScriptName = self.pcssRunner.internalConfig["training_annotation_node_script"]
        shellScriptName = self.pcssRunner.internalConfig["training_annotation_shell_script"]
        outputFileName = self.pcssRunner.internalConfig["training_annotation_stdout_file"]
        self.makeFullAnnotationSgeScript(nodeScriptName, shellScriptName, outputFileName)       

    def makeFullAnnotationSgeScript(self, nodeScriptName, shellScriptName, outputFileName):
        script = self.csg.makeClusterHeaderCommands(self.seqBatchCount)
        script += self.csg.makeBaseAnnotationSgeScript(self.getSequenceTaskListString(), self.getSeqBatchDirectory(), nodeScriptName, outputFileName)
        scriptOutputFile = self.pcssRunner.pdh.getFullOutputFile(shellScriptName)
        scriptFh = open(scriptOutputFile, 'w')
        scriptFh.write(script)
        scriptFh.close()

class ClusterScriptGenerator:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner

    def makeClusterHeaderCommands(self, taskCount):

        script = """
#!/bin/tcsh
#$ -S /bin/tcsh                    
#$ -o clusterOutput.txt            
#$ -e clusterError.txt             
#$ -cwd                            
#$ -r y                            
#$ -j y                            
#$ -l mem_free=1G                  
#$ -l arch=linux-x64               
#$ -l netapp=1G,scratch=1G         
#$ -l h_rt=24:00:00                
#$ -t 1-%(taskCount)s

""" %locals()
        return script

    def makeBaseAnnotationSgeScript(self, taskList, seqBatchDir, commandName, outputFileName):

        pcssBaseDirectory = self.pcssRunner.pcssConfig["pcss_directory"]
                
        parameterFileName = self.pcssRunner.internalConfig["seq_batch_parameter_file_name"]
        nodeHomeDirectory = self.pcssRunner.internalConfig["cluster_pipeline_directory"]
        
        taskListString = " ".join(taskList)

        script = """

set PCSS_BASE_DIRECTORY="%(pcssBaseDirectory)s"

set MODEL_OUTPUT_FILE_NAME="%(outputFileName)s"

set tasks=( %(taskListString)s )
set input=$tasks[$SGE_TASK_ID]

set NODE_HOME_DIR="%(nodeHomeDirectory)s/$input"
mkdir -p $NODE_HOME_DIR
cd $NODE_HOME_DIR

date
hostname
pwd

set PARAMETER_FILE_NAME="%(seqBatchDir)s/$input/%(parameterFileName)s"

setenv PYTHONPATH $PCSS_BASE_DIRECTORY/lib
python $PCSS_BASE_DIRECTORY/bin/clusterExe/%(commandName)s $PARAMETER_FILE_NAME > & $MODEL_OUTPUT_FILE_NAME

cp $MODEL_OUTPUT_FILE_NAME "%(seqBatchDir)s/$input/"

rm -r $NODE_HOME_DIR/
""" %locals()

        return script

    def makeBaseBenchmarkSgeScript(self):

        pcssBaseDirectory = self.pcssRunner.pcssConfig["pcss_directory"]
                
        parameterFileName = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["benchmark_parameter_file_name"])
        modelOutputFileName = self.pcssRunner.internalConfig["training_benchmark_stdout_file"]
        nodeHomeDirectory = self.pcssRunner.internalConfig["cluster_pipeline_directory"]
        commandName = self.pcssRunner.internalConfig["training_benchmark_node_script"]

        headNodeOutputDirectory = self.pcssRunner.pdh.getFullOutputFile("")
        
        script = """

set PCSS_BASE_DIRECTORY="%(pcssBaseDirectory)s"

set MODEL_OUTPUT_FILE_NAME="%(modelOutputFileName)s"

set tasks=( 1 )
set input=$tasks[$SGE_TASK_ID]

set NODE_HOME_DIR="%(nodeHomeDirectory)s/$input"
mkdir -p $NODE_HOME_DIR
cd $NODE_HOME_DIR

date
hostname
pwd

set PARAMETER_FILE_NAME="%(parameterFileName)s"

setenv PYTHONPATH $PCSS_BASE_DIRECTORY/lib
python $PCSS_BASE_DIRECTORY/bin/clusterExe/%(commandName)s $PARAMETER_FILE_NAME > & $MODEL_OUTPUT_FILE_NAME

cp $MODEL_OUTPUT_FILE_NAME "%(headNodeOutputDirectory)s/"

rm -r $NODE_HOME_DIR/
""" %locals()

        return script


class ClusterBenchmarker:

    def __init__(self, pcssRunner):
        self.pdh = pcssRunner.pdh
        self.pcssRunner = pcssRunner
        self.csg = ClusterScriptGenerator(pcssRunner)

    def prepareTrainingBenchmarkRun(self):
        pcssCopy = copy.deepcopy(self.pcssRunner.pcssConfig)
        print "running cluster with output file %s" % self.pcssRunner.pdh.getFullOutputFile("")
        pcssCopy["using_web_server"] = False
        inputAnnotationFileName = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["annotation_output_file"])
        if (not os.path.exists(inputAnnotationFileName)):
            msg = "Did not find input annotation file name %s\n" % inputAnnotationFileName
            msg += "Please make sure this is file is in the run directory for this training benchmark run"
            raise pcssErrors.PcssGlobalException(msg)
            
        pcssCopy["input_annotation_file_name"] = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["annotation_output_file"])
        pcssCopy.filename = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["benchmark_parameter_file_name"])
        pcssCopy.write()
        
    def makeFullTrainingBenchmarkScript(self):
        script = self.csg.makeClusterHeaderCommands(1)
        script += self.csg.makeBaseBenchmarkSgeScript()
        scriptOutputFile = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["training_benchmark_shell_script"])
        scriptFh = open(scriptOutputFile, 'w')
        scriptFh.write(script)
        scriptFh.close()

