import sys
import os
from Bio import SeqIO
import itertools
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
        self.pcssRunner.setJobDirectory(os.path.join(self.pcssRunner.pcssConfig["run_directory"], "developClusterJob"))
        self.subDirectories = []

    def getFastaGroupList(self, fastaFile):
        fh = open(fastaFile, 'r')
        fastaList = list(SeqIO.FastaIO.FastaIterator(fh))
        seqBatchSize = int(self.pcssRunner.internalConfig["seq_batch_size"])
        print "seq batch size: %s" % seqBatchSize
        seqGroupList = [fastaList[i:i+seqBatchSize] for i in range(0, len(fastaList), seqBatchSize)]
        fh.close()
        self.seqBatchCount = len(seqGroupList)
        return seqGroupList

    def getSequenceTaskList(self):
        taskList = []
        for i in range(self.seqBatchCount):
            taskList.append(str(i))
        return " ".join(taskList)
                

    #when we get here, we have fasta file and unannotated sequences. 
    def divideSeqsFromFasta(self, fastaFile):
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
        
    def makeBaseSgeScript(self):

        taskList = self.getSequenceTaskList()

        pcssBaseDirectory = self.pcssRunner.pcssConfig["pcss_directory"]
        topLevelSeqBatchDir = self.getSeqBatchDirectory()
                
        parameterFileName = self.pcssRunner.internalConfig["seq_batch_parameter_file_name"]
        modelOutputFileName = self.pcssRunner.internalConfig["cluster_stdout_file"]
        nodeHomeDirectory = self.pcssRunner.internalConfig["cluster_pipeline_directory"]
        modelPipelineScriptName = self.pcssRunner.internalConfig["model_pipeline_script_name"]
        
        taskListString = " ".join(taskList)

        script = """

set PCSS_BASE_DIRECTORY="%(pcssBaseDirectory)s"

set MODEL_OUTPUT_FILE_NAME="%(modelOutputFileName)s"

set tasks=( %(taskListString)s )
set input=$tasks[$SGE_TASK_ID]

set NODE_HOME_DIR="%(nodeHomeDirectory)s/$input"
mkdir -p $NODE_HOME_DIR
cd $NODE_HOME_DIR

date
hostname
pwd

set PARAMETER_FILE_NAME="%(topLevelSeqBatchDir)s/$input/%(parameterFileName)s"

setenv PYTHONPATH $PCSS_BASE_DIRECTORY/lib
python $PCSS_BASE_DIRECTORY/bin/%(modelPipelineScriptName)s --parameterFileName $PARAMETER_FILE_NAME  > & $MODEL_OUTPUT_FILE_NAME

rm -r $NODE_HOME_DIR/
""" %locals()

        scriptOutputFile = self.pcssRunner.pdh.getFullOutputFile("develop_cluster_script.sh")
        scriptFh = open(scriptOutputFile, 'w')
        scriptFh.write(script)
        scriptFh.close()

        return script

