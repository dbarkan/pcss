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
        self.pcssRunner = pcssRunner
        self.pdh = pcssRunner.pdh
        self.subDirectories = []

    def getFastaGroupList(self, fastaFile):
        fh = open(fastaFile, 'r')
        fastaList = sorted(list(SeqIO.FastaIO.FastaIterator(fh)), key=lambda x: x.id)
        seqBatchSize = int(self.pcssRunner.internalConfig["seq_batch_size"])

        print "seq batch size: %s" % seqBatchSize
        seqGroupList = [fastaList[i:i+seqBatchSize] for i in range(0, len(fastaList), seqBatchSize)]
        fh.close()
        self.seqBatchCount = len(seqGroupList)
        return seqGroupList

    def getSeqBatchCount(self):
        return self.seqBatchCount

    def getSequenceTaskListString(self):
        taskList = []
        for i in range(self.seqBatchCount):
            taskList.append(str(i))
        return " ".join(taskList)
 
    def mergeSvmApplicationResults(self):
        
        fastaFile = self.pcssRunner.pcssConfig['fasta_file']
        seqGroupList = self.getFastaGroupList(fastaFile)
        print "merging application results"
        seqBatchDirectory = self.pdh.getSeqBatchDirectory()
        allProteins = []
        for (i, nextGroup) in enumerate(seqGroupList):
            subDirName = self.pdh.getSeqBatchSubDirectoryName(i)
            self.seqBatchErrorExists(subDirName)
            subOutputFile = os.path.join(subDirName, self.pcssRunner.internalConfig["annotation_output_file"])
            if (not os.path.exists(subOutputFile)):
                raise pcssErrors.PcssGlobalException("Seq batch error: did not get annotation output file in directory %s" % subDirName)
            reader = pcssIO.AnnotationFileReader(self.pcssRunner)
            reader.readAnnotationFile(subOutputFile)
            proteins = reader.getProteins()
            allProteins += proteins
        
        afw = pcssIO.AnnotationFileWriter(self.pcssRunner)
        afw.writeAllOutput(allProteins)

    def seqBatchErrorExists(self, subDirName):
        if (os.path.exists(os.path.join(subDirName, self.pcssRunner.internalConfig["pcss_error_output_file"]))):
            errorInfo = pcssErrors.ErrorInfo(os.path.join(subDirName, self.pcssRunner.internalConfig["pcss_error_output_file"]))
            raise pcssErrors.PcssGlobalException("Got pcss seq batch error %s\nin directory %s" % (errorInfo.msg, subDirName))

        if (os.path.exists(os.path.join(subDirName, self.pcssRunner.internalConfig["internal_error_output_file"]))):
            errorInfo = pcssErrors.ErrorInfo(os.path.join(subDirName, self.pcssRunner.internalConfig["internal_error_output_file"]))
            raise InternalException("Got internal seq batch error %s\nin directory %s" % (errorInfo.msg, subDirName))

    def mergeTrainingAnnotationResults(self):
        
        fastaFile = self.pcssRunner.pcssConfig['fasta_file']
        seqGroupList = self.getFastaGroupList(fastaFile)

        seqBatchDirectory = self.pdh.getSeqBatchDirectory()
        allProteins = []
        for (i, nextGroup) in enumerate(seqGroupList):
            subDirName = self.pdh.getSeqBatchSubDirectoryName(i)
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

        seqBatchDirectory = self.pdh.getSeqBatchDirectory()
        
        for (i, nextGroup) in enumerate(seqGroupList):
            self.createSequenceSubDirectory(i)

            self.makeSeqBatchFile(i, nextGroup)
        
    def createSequenceSubDirectory(self, i):
        subDirectoryName = self.pdh.getSeqBatchSubDirectoryName(i)
        if (not os.path.exists(subDirectoryName)):
            os.mkdir(subDirectoryName)
        
    def makeSeqBatchFile(self, i, seqBatch):
        fileName = self.pdh.getSeqBatchFastaFileName(i)
        fh = open(fileName, 'w')
        for seqRecord in seqBatch:
            fh.write(">%s\n" % seqRecord.id)
            fh.write("%s\n" % str(seqRecord.seq))
        fh.close()

    def makeRunDisopredSgeScript(self):
        nodeScriptName = self.pcssRunner.internalConfig["disopred_standalone_node_script"]
        shellScriptName = self.pcssRunner.internalConfig["disopred_standalone_shell_script"]
        outputFileName = self.pcssRunner.internalConfig["disopred_standalone_stdout_file"]
        self.makeFullAnnotationSgeScript(nodeScriptName, shellScriptName, outputFileName)
        
    def makeFullTrainingAnnotationSgeScript(self):
        nodeScriptName = self.pcssRunner.internalConfig["training_annotation_node_script"]
        shellScriptName = self.pcssRunner.internalConfig["training_annotation_shell_script"]
        outputFileName = self.pcssRunner.internalConfig["training_annotation_stdout_file"]
        self.makeFullAnnotationSgeScript(nodeScriptName, shellScriptName, outputFileName)       

class ConfigFileGenerator:

    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.pdh = pcssRunner.pdh

    def setSeqDivider(self, seqDivider):
        self.seqDivider = seqDivider

class SvmApplicationConfigFileGenerator(ConfigFileGenerator):
    def generateConfigFiles(self):
        for i in range(self.seqDivider.getSeqBatchCount()):
            pcssConfig = self.pdh.getClusterNodeConfig(i)

            subDirectoryName = self.pdh.getSeqBatchSubDirectoryName(i)
            pcssConfig.filename = os.path.join(subDirectoryName, self.pcssRunner.internalConfig["seq_batch_node_config_file"])
            pcssConfig.write()
    
class ClusterScriptGenerator:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.pdh = self.pcssRunner.pdh

    def writeScript(self, script, shellScriptName):
        scriptOutputFile = self.pcssRunner.pdh.getFullOutputFile(shellScriptName)
        scriptFh = open(scriptOutputFile, 'w')
        scriptFh.write(script)
        scriptFh.close()
        
    def makeClusterHeaderCommands(self, taskCount):

        script = """
#!/bin/tcsh
#$ -S /bin/tcsh                    
#$ -o sgeOutput.txt            
#$ -e sgeError.txt             
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

    def makeBaseAnnotationSgeScript(self, taskListString, commandName, outputFileName):

        pcssClusterBaseDirectory = self.pdh.getPcssClusterBaseDirectory()
        seqBatchDir = self.pdh.getClusterSeqBatchDirectory()
        configFileName = self.pcssRunner.internalConfig["seq_batch_node_config_file"]
        nodeHomeDirectory = self.pcssRunner.internalConfig["cluster_pipeline_directory"]
        
        script = """

set PCSS_BASE_DIRECTORY="%(pcssClusterBaseDirectory)s"

set MODEL_OUTPUT_FILE_NAME="%(outputFileName)s"

set tasks=( %(taskListString)s )
set input=$tasks[$SGE_TASK_ID]

set NODE_HOME_DIR="%(nodeHomeDirectory)s/$input"
mkdir -p $NODE_HOME_DIR
cd $NODE_HOME_DIR

date
hostname
pwd

set CONFIG_FILE_NAME="%(seqBatchDir)s/$input/%(configFileName)s"

setenv PYTHONPATH $PCSS_BASE_DIRECTORY/lib
python $PCSS_BASE_DIRECTORY/bin/clusterExe/%(commandName)s $CONFIG_FILE_NAME > & $MODEL_OUTPUT_FILE_NAME

cp $MODEL_OUTPUT_FILE_NAME "%(seqBatchDir)s/$input/"

rm -r $NODE_HOME_DIR/
""" %locals()

        return script

    def makeBaseBenchmarkSgeScript(self):

        pcssBaseDirectory = self.pcssRunner.pcssConfig["pcss_directory"]
                
        configFileName = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["training_benchmark_config_file"])
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

set CONFIG_FILE_NAME="%(configFileName)s"

setenv PYTHONPATH $PCSS_BASE_DIRECTORY/lib
python $PCSS_BASE_DIRECTORY/bin/clusterExe/%(commandName)s $CONFIG_FILE_NAME > & $MODEL_OUTPUT_FILE_NAME

cp $MODEL_OUTPUT_FILE_NAME "%(headNodeOutputDirectory)s/"

rm -r $NODE_HOME_DIR/
""" %locals()

        return script

class SvmApplicationClusterScriptGenerator(ClusterScriptGenerator):

    def makeBaseSvmApplicationSgeScript(self):
        taskListString = self.seqDivider.getSequenceTaskListString() 
        nodeScriptName = self.pcssRunner.internalConfig["svm_application_node_script"]
        outputFileName = self.pcssRunner.internalConfig["svm_application_stdout_file"]
        return self.makeBaseAnnotationSgeScript(taskListString, nodeScriptName, outputFileName)

    def writeFullSvmApplicationSgeScript(self):
        script = self.makeClusterHeaderCommands(self.seqDivider.getSeqBatchCount())
        script += self.makeBaseSvmApplicationSgeScript()
        shellScriptName = self.pcssRunner.internalConfig["svm_application_shell_script"]
        self.writeScript(script, shellScriptName)

    def setSeqDivider(self, seqDivider):
        self.seqDivider = seqDivider

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
        pcssCopy.filename = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["training_benchmark_config_file"])
        pcssCopy.write()
        
    def makeFullTrainingBenchmarkScript(self):
        script = self.csg.makeClusterHeaderCommands(1)
        script += self.csg.makeBaseBenchmarkSgeScript()
        scriptOutputFile = self.pcssRunner.pdh.getFullOutputFile(self.pcssRunner.internalConfig["training_benchmark_shell_script"])
        scriptFh = open(scriptOutputFile, 'w')
        scriptFh.write(script)
        scriptFh.close()


class TrainingAnnotationClusterScriptGenerator(ClusterScriptGenerator):

    def makeBaseTrainingAnnotationSgeScript(self, taskListString):
        nodeScriptName = self.pcssRunner.internalConfig["training_annotation_node_script"]
        outputFileName = self.pcssRunner.internalConfig["training_annotation_stdout_file"]
        return self.makeBaseAnnotationSgeScript(taskListString, nodeScriptName, outputFileName)

    def writeFullTrainingAnnotationSgeScript(self):
        script = self.makeClusterHeaderCommands(self.seqDivider.getSeqBatchCount())
        script += self.makeBaseTrainingAnnotationSgeScript(self.seqDivider.getSequenceTaskListString())
        shellScriptName = self.pcssRunner.internalConfig["training_annotation_shell_script"]
        self.writeScript(script, shellScriptName)

    def setSeqDivider(self, seqDivider):
        self.seqDivider = seqDivider

