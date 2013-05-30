import sys
from validate import ValidateError
import os
import re
from Bio import SeqIO
import collections
import myCollections
import pcssIO
import pcssCluster
import tempfile
import configobj
import copy
import subprocess
import time
import logging
import pcssSvm
import pcssErrors
from configobj import flatten_errors
from validate import Validator
import pcssModels
import pcssFeatures
import pcssFeatureHandlers
import shutil
import traceback
from Bio import PDB
log = logging.getLogger("pcssTools")

def getOneLetterFromBioResidue(residueObject):
    return PDB.Polypeptide.three_to_one(residueObject)

def getProteinErrorCodePrefix():
    return "protein_error_"

def getPeptideErrorCodePrefix():
    return "peptide_error_"

def isPeptideErrorValue(value):
    return value.startswith(getPeptideErrorCodePrefix())

def getAllPeptides(proteins, ignoreProteinErrors):
    allPeptides = []
    for protein in proteins:
        if (ignoreProteinErrors or not protein.hasErrors()):
            allPeptides = allPeptides + protein.peptides.values()
    return allPeptides

def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def fileExists(fileName):
    if (not os.path.exists(fileName)):
        raise ValidateError("FileName %s does not exist" % fileName)
    return fileName

class PcssRunner:

    """Class for running PCSS algorithms. Handles directory structure and provides handle to config object, among others"""

    def __init__(self, pcssConfig):
        """Create new PTP runner. Initialization creates PCSS Directory Handler"""
        self.pcssConfig = pcssConfig
        if (pcssConfig.configspec is not None):
            self.validateConfig(pcssConfig)

        internalConfigFile = os.path.join(self.pcssConfig["user_pcss_directory"], "data", "config", "internalConfigFile.txt")
        internalConfigSpec = os.path.join(self.pcssConfig["user_pcss_directory"], "data", "config", "internalConfigSpec.txt")

        self.internalConfig = configobj.ConfigObj(internalConfigFile, configspec=internalConfigSpec)
        self.validateConfig(self.internalConfig)

        self.pdh = self.createDirectoryHandler(pcssConfig, self.internalConfig)
        self.pdh.createOutputDirectory()

        self.parser = PDB.PDBParser(QUIET=True)
        self.readFileAttributes()
        self.peptideLength = None
        
        logging.basicConfig(filename=self.pdh.getFullOutputFile("%s.log" % self.getRunName()), level=logging.DEBUG,
                            filemode="w", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

        self.initSubclass()

    def createDirectoryHandler(self, pcssConfig, internalConfig):
        return PcssDirectoryHandler(pcssConfig, internalConfig)

    def initSubclass(self):
        return

    def execute(self):
        try:
            log.info("Beginning pcss run")
            self.checkForErrors()
            self.executePipeline()
            log.info("Finished pcss run")
        except pcssErrors.PcssGlobalException, pge:
            self.handlePcssGlobalException(pge)

        except pcssErrors.ErrorExistsException, eee:
            self.handleErrorExistsException(eee)

        except pcssErrors.InternalException, e:
            self.handleInternalException(e)
        
        except Exception, e:
            self.handlePythonException(e)


    def validateConfig(self, config):
        validator = Validator({'file': fileExists})
        results = config.validate(validator, preserve_errors=True)
        if (results != True):
            self.handleConfigError(results)

    def readFileAttributes(self):
        return

            
    def checkForErrors(self):
        if (os.path.exists(self.pdh.getPcssErrorFile())):
            raise pcssErrors.ErrorExistsException("Previous error exists", self.pdh.getPcssErrorFile())
        if (os.path.exists(self.pdh.getInternalErrorFile())):
            raise pcssErrors.ErrorExistsException("Previous error exists", self.pdh.getInternalErrorFile())

    def validatePeptideCodeStatus(self, status, peptideCode):
        try:
            return self.validatePeptideStatus(status)
        except pcssErrors.PcssGlobalException, e:
            raise pcssErrors.PcssGlobalException("%s (peptide code: %s)" % (e.msg, peptideCode))

    def validatePeptideStatus(self, status):
                                                
        status = status.lower()
        if (not (status == self.getPositiveKeyword() or status == self.getNegativeKeyword() or status == self.getApplicationKeyword())):
            raise pcssErrors.PcssGlobalException("Peptide status %s not valid status (needs to be %s, %s, %s" % (status, self.getPositiveKeyword(), 
                                                                                                                 self.getNegativeKeyword(), self.getApplicationKeyword()))
        return status

    def validatePeptideTrainingStatus(self, status):
                                                
        status = status.lower()
        if (not (status == self.getPositiveKeyword() or status == self.getNegativeKeyword())):
            raise pcssErrors.PcssGlobalException("Peptide status %s not valid status (needs to be %s, %s)" % (status, self.getPositiveKeyword(), 
                                                                                                              self.getNegativeKeyword()))
        return status

    def setPeptideLength(self, peptideLength):
        self.peptideLength = peptideLength

    def getPeptideLength(self):
        if (self.peptideLength is None):
            raise pcssErrors.PcssGlobalException("Error: Peptide Length was never set; please make sure it is set when reading input")
        return self.peptideLength

    def getBioParser(self):
        return self.parser

    def getSvmFeatureOrder(self):
        featureOrderList = self.internalConfig["feature_order"]
        return featureOrderList

    def writeErrorFile(self, errorType, message, fileName):
        errorFh = open(fileName, 'w')
        errorFh.write(errorType + "\n")
        errorFh.write(message + "\n")
        tb = traceback.format_exc()
        errorFh.write(tb + "\n")
        print message
        print tb
        log.error(message)
        log.error(tb)
        errorFh.close()
    def handlePcssGlobalException(self, pge):
        errorFile = self.pdh.getPcssErrorFile()
        self.writeErrorFile(self.internalConfig["keyword_pcss_error"], "PCSS ERROR: " + pge.msg, errorFile)

    def handleErrorExistsException(self, eee):
        errorInfo = pcssErrors.ErrorInfo(eee.fileName)
        msg = "EXISTING ERROR: A previous run in the same job directory finished with the following error\n"
        msg += errorInfo.msg
        print msg
        log.error(msg)

    def handleInternalErrorFile(self, e):
        errorFile = self.pdh.getInternalErrorFile()
        self.writeErrorFile(self.internalConfig["keyword_internal_error"], "INTERNAL ERROR: " + e.msg, errorFile)

    def handlePythonException(self, e):
        errorFile = self.pdh.getInternalErrorFile()
        self.writeErrorFile(self.internalConfig["keyword_internal_error"], "INTERNAL ERROR: " + str(e), errorFile)

    def getErrorFileInfo(self):
        if (os.path.exists(self.pdh.getInternalErrorFile())):
            errorInfo = pcssErrors.ErrorInfo(self.pdh.getInternalErrorFile())
            return errorInfo

        elif (os.path.exists(self.pdh.getPcssErrorFile())):
            errorInfo = pcssErrors.ErrorInfo(self.pdh.getPcssErrorFile())
            return errorInfo
        else:
            return None


    def getPositiveKeyword(self):
        return self.internalConfig["keyword_positive_status"]

    def getNegativeKeyword(self):
        return self.internalConfig["keyword_negative_status"]

    def getApplicationKeyword(self):
        return self.internalConfig["keyword_application_status"]

    def readAnnotationFile(self):
        self.reader = pcssIO.AnnotationFileReader(self)
        self.reader.readAnnotationFile(self.pcssConfig["input_annotation_file_name"])
        self.proteins = self.reader.getProteins()
        
    def getModelUrl(self, modelId):
        modelUrlInitial = self.internalConfig["model_url"]
        wildCard = self.internalConfig["wild_card"]
        modelUrlFinal = modelUrlInitial.replace(wildCard, modelId)
        return modelUrlFinal

    def getRunName(self):
        return self.pdh.getRunName()

    def readProteins(self):

        peptideImporterType = self.pcssConfig["peptide_importer_type"]
        peptideImporter = None
        if (peptideImporterType == "scan"):
            peptideImporter = pcssIO.ScanPeptideImporter(self)
        elif(peptideImporterType == "defined"):
            
            peptideImporter = pcssIO.DefinedPeptideImporter(self)
        else:
            peptideImporter = pcssIO.FullProteinImporter(self)
        self.proteins = peptideImporter.readInputFile(self.pcssConfig['fasta_file'])

    def handleConfigError(self, results):
        msg = "CONFIGURATION ERROR\n"
        for (section_list, key, _) in flatten_errors(self.pcssConfig, results):
            if key is not None:
                msg +=  'The "%s" key in the section "%s" failed validation\n' % (key, ', '.join(section_list))
            else:
                msg += 'The following section was missing:%s ' % ', '.join(section_list)
        print msg
        raise pcssErrors.PcssGlobalException(msg)
                
    def addPeptideFeatures(self):
        
        modelColumns = pcssModels.PcssModelTableColumns(self.internalConfig['model_table_column_file'])
        modelTable = pcssModels.PcssModelTable(self, modelColumns)

        disopredFileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.pdh)
        disopredReader = pcssFeatureHandlers.DisopredReader(disopredFileHandler)
        disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(disopredFileHandler)
        
        psipredFileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.pdh)
        psipredReader = pcssFeatureHandlers.PsipredReader(psipredFileHandler)
        psipredRunner = pcssFeatureHandlers.SequenceFeatureRunner(psipredFileHandler)

        for protein in self.proteins:
            protein.processDisopred(disopredReader, disopredRunner)
            protein.processPsipred(psipredReader, psipredRunner)
            protein.addModels(modelTable)        
            protein.processDssp()
            
        #sum1 = summary.summarize(muppy.get_objects())
        #summary.print_(sum1)
    def writeOutput(self):

        afw = pcssIO.AnnotationFileWriter(self)
        afw.writeAllOutput(self.proteins)

class PrepareDisopredClusterRunner(PcssRunner):
    def executePipeline(self):
        seqDivider = pcssCluster.SeqDivider(self)
        
        seqDivider.divideSeqsFromFasta()

        seqDivider.makeRunDisopredSgeScript()

class DisopredStandaloneRunner(PcssRunner):
    
    def executePipeline(self):

        self.readProteins()

        self.runDisopred()

    def runDisopred(self):

        disopredFileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.pdh)
        disopredReader = pcssFeatureHandlers.DisopredReader(disopredFileHandler)
        disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(disopredFileHandler)

        for protein in self.proteins:
            protein.processDisopred(disopredReader, disopredRunner)

class ModelRunner(PcssRunner):
    def initSubclass(self):
        self.modelHandler = PcssModelHandler(self.pcssConfig, self.pdh)

class PrepareClusterRunner(PcssRunner):
    def updateInputFileConfig(self):
        self.pcssConfig["fasta_file"] = self.pdh.getFullOutputFile(self.internalConfig["server_input_fasta_file_name"])
        self.pcssConfig["rules_file"] = self.pdh.getFullOutputFile(self.internalConfig["server_input_rules_file_name"])

    def prepareDirectories(self, cfg):
        seqDivider = pcssCluster.SeqDivider(self)
        
        seqDivider.divideSeqsFromFasta()
        
        cfg.setSeqDivider(seqDivider)

        cfg.generateConfigFiles()
        
        self.seqDivider = seqDivider

    def getTaskCount(self):
        return self.seqDivider.getSeqBatchCount()

    def setClusterShellScript(self, script):
        self.clusterShellScript = script
        
    def getClusterShellScript(self):
        return self.clusterShellScript


class PrepareSvmApplicationClusterRunner(PrepareClusterRunner):
    def executePipeline(self):

        cfg = pcssCluster.SvmApplicationConfigFileGenerator(self)
        
        self.prepareDirectories(cfg)
        
        csg = pcssCluster.SvmApplicationClusterScriptGenerator(self)

        csg.setSeqDivider(self.seqDivider)

        csg.writeFullSvmApplicationSgeScript()

class PrepareTrainingBenchmarkClusterRunner(PcssRunner):

    def executePipeline(self):
        tcb =  pcssCluster.ClusterBenchmarker(self)
            
        tcb.prepareTrainingBenchmarkRun()
            
        tcb.makeFullTrainingBenchmarkScript()
        
class PrepareTrainingAnnotationClusterRunner(PrepareClusterRunner):
    def executePipeline(self):
        cfg = pcssCluster.SvmApplicationConfigFileGenerator(self)
        
        self.prepareDirectories(cfg)
        
        csg = pcssCluster.TrainingAnnotationClusterScriptGenerator(self)

        csg.setSeqDivider(self.seqDivider)

        csg.writeFullTrainingAnnotationSgeScript()
            
class FinalizeApplicationClusterRunner(PcssRunner):
    def executePipeline(self):
        seqDivider = pcssCluster.SeqDivider(self)

        seqDivider.mergeSvmApplicationResults()

    def readFileAttributes(self):
        fileName = self.internalConfig["svm_application_cluster_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)


class FinalizeApplicationServerRunner(FinalizeApplicationClusterRunner):
    def createDirectoryHandler(self, pcssConfig, internalConfig):
        return PcssServerDirectoryHandler(pcssConfig, internalConfig)

    def executePipeline(self):
        self.updateInputFileConfig()
        seqDivider = pcssCluster.SeqDivider(self)
        seqDivider.mergeSvmApplicationResults()
    def updateInputFileConfig(self):
        self.pcssConfig["fasta_file"] = self.pdh.getFullOutputFile(self.internalConfig["server_input_fasta_file_name"])
        self.pcssConfig["rules_file"] = self.pdh.getFullOutputFile(self.internalConfig["server_input_rules_file_name"])


class PrepareTrainingAnnotationServerRunner(PrepareTrainingAnnotationClusterRunner):
    def executePipeline(self):
        self.updateInputFileConfig()

        cfg = pcssCluster.SvmApplicationConfigFileGenerator(self)

        self.prepareDirectories(cfg)

        csg = pcssCluster.TrainingAnnotationClusterScriptGenerator(self)

        csg.setSeqDivider(self.seqDivider)

        self.setClusterShellScript(csg.makeBaseSvmApplicationSgeScript())

    def createDirectoryHandler(self, pcssConfig, internalConfig):
        return PcssServerDirectoryHandler(pcssConfig, internalConfig)

      
class PrepareSvmApplicationServerRunner(PrepareSvmApplicationClusterRunner):
    def executePipeline(self):
        self.updateInputFileConfig()

        cfg = pcssCluster.SvmApplicationConfigFileGenerator(self)

        self.prepareDirectories(cfg)

        csg = pcssCluster.SvmApplicationClusterScriptGenerator(self)

        csg.setSeqDivider(self.seqDivider)

        self.setClusterShellScript(csg.makeBaseSvmApplicationSgeScript())
  
    def createDirectoryHandler(self, pcssConfig, internalConfig):
        return PcssServerDirectoryHandler(pcssConfig, internalConfig)

class SvmApplicationRunner(ModelRunner):
    def runSvm(self):
        self.appSvm = pcssSvm.ApplicationSvm(self)
        self.appSvm.setProteins(self.proteins)
        self.appSvm.writeClassificationFile()
        self.appSvm.classifySvm()
        self.appSvm.readResultFile()
        self.appSvm.addScoresToPeptides()

    def readFileAttributes(self):
        fileName = self.internalConfig["svm_application_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)


class SvmApplicationInputRunner(SvmApplicationRunner):

    def executePipeline(self):

        self.readSvmApplicationInput()

        self.runSvm()
        
        self.writeOutput()
        
    def readSvmApplicationInput(self):
        self.readAnnotationFile()
        
        self.reader.readProteinSequences(self.pcssConfig["fasta_file"])

class SvmApplicationFeatureRunner(SvmApplicationRunner):

    def executePipeline(self):

        self.readProteins()

        self.addPeptideFeatures()

        self.runSvm()
        
        self.writeOutput()
        
class AnnotationRunner(ModelRunner):

    def executePipeline(self):
        self.readProteins()

        self.addPeptideFeatures()
        
        self.writeOutput()

    def readFileAttributes(self):

        fileName = self.internalConfig["annotation_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)

class TrainingAnnotationRunner(AnnotationRunner):
    def readFileAttributes(self):

        fileName = self.internalConfig["training_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)
    

class CompleteSvmRunner(PcssRunner):
    def executePipeline(self):
        self.readAnnotationFile()
        
        self.createSvmModel()

    def createSvmModel(self):
        generator = pcssSvm.CompleteSvmGenerator(self)

        generator.createSvmModel(getAllPeptides(self.reader.getProteins(), False))

    def readFileAttributes(self):

        fileName = self.internalConfig["training_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)
    
class LeaveOneOutBenchmarkRunner(PcssRunner):
    def executePipeline(self):
        self.readAnnotationFile()
        
        self.benchmark()

    def benchmark(self):
        benchmarker = pcssSvm.LeaveOneOutBenchmarker(self)
        peptideSet = getAllPeptides(self.reader.getProteins(), False)
        for i in range(len(peptideSet)):
            benchmarker.createTrainingAndTestSets(peptideSet)
            benchmarker.trainAndApplyModel() 
            benchmarker.readBenchmarkResults()

        benchmarker.processAllResults()

    def readFileAttributes(self):

        fileName = self.internalConfig["training_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)
        
class TrainingBenchmarkRunner(PcssRunner):
    def executePipeline(self):

        self.readAnnotationFile()
        
        self.benchmark()
              
    def benchmark(self):
        
        benchmarker = pcssSvm.SvmBenchmarker(self)
 
        for i in range(int(self.pcssConfig["training_iterations"])):
            benchmarker.createTrainingAndTestSets(getAllPeptides(self.proteins, False))  #parition pepitdes, write training and test set files

            benchmarker.trainAndApplyModel() #call svm command line

            benchmarker.readBenchmarkResults() #read result file

        benchmarker.processAllResults()

    def readFileAttributes(self):

        fileName = self.internalConfig["training_attribute_file"]
        self.pfa = pcssIO.PcssFileAttributes(fileName)
    

class PcssModelHandler:

    """Class for managing model PDB files, retrieving them from file servers as necessary"""

    def __init__(self, pcssConfig, pdh):
        self.pcssConfig = pcssConfig
        self.pdh = pdh
        self.loadRunInfo()

    def loadRunInfo(self):
        """Load PcssModelRunInfo object which provide data specific to differen modpipe runs"""
        
        self._runInfo = {}
        runInfoFile  = self.pdh.internalConfig["model_run_info"]
        runInfoConfig = configobj.ConfigObj(runInfoFile)
        
        for runName, values in runInfoConfig.iteritems():
            self._runInfo[runName] = PcssModelRunInfo(runName, values, self.pdh)
                         
    def getLocalModelFileName(self, pcssModel):
        """Get model file name as stored where PCSS expects; if not present, search on modbase file server"""
        fullModelFileName = self.pdh.getFullModelFile(pcssModel)
        log.debug("Checking for existence model file %s" %  fullModelFileName)
        if (not os.path.exists(fullModelFileName)):
            log.debug("Model file %s doesn't exist; retrieving from source directory" % fullModelFileName)
            self.retrieveModelFile(pcssModel)
        return fullModelFileName

    def retrieveModelFile(self, pcssModel):
        """Copy model file from modbase file server to local pcss directory"""
        modelRunInfo = self._runInfo[pcssModel.getRunName()]
        sourcePath = modelRunInfo.getSourcePath(pcssModel)
        sourceModelFileZip = self.makeSourceModelZipFileName(pcssModel)
        fullSourceModelFileZip = os.path.join(sourcePath, sourceModelFileZip)
        if (not(os.path.exists(fullSourceModelFileZip))):
            raise pcssErrors.NoSourceModelException("Did not find model file in source directory (searched for %s" % fullSourceModelFileZip)

        self.pdh.copyFile(sourcePath, sourceModelFileZip, self.pdh.getStructureDirectory())
        self.pdh.unzipFile(os.path.join(self.pdh.getStructureDirectory(), sourceModelFileZip))

    def makeSourceModelZipFileName(self, model):
        return "%s.pdb.gz" % model.getId()

class PcssModelRunInfo:

    """Simple class tracking data specific to modpipe runs; for example, each run has its own directory where models are written"""

    def __init__(self, name, values, pdh):
        self.path = values['path']
        self.offset = values['offset']
        self.style = self.createModelStyle(name, values, pdh)
        self.pdh = pdh

    def createModelStyle(self, runName, values, pdh):
        if (values['style'] == 'NewModelStyle'):
            return NewModelStyle(pdh)
        elif (values['style'] == 'OldModelStyle'):
            return OldModelStyle(pdh)
        else:
            raise pcssErrors.PcssGlobalException("Model run info file has run %s with invalid model style %s; "
                                                 "please change to either 'NewModelStyle' or 'OldModelStyle'" %
                                                 (runName, values['style']))
    def getSourcePath(self, pcssModel):
        return self.style.getSourcePath(self.path, pcssModel)

class ModelStyle:


    """Base class providing information specific to model run style.

    Around 2006 the way modbase models were stored along with their residue numbering changed. This 
    class manages the differences"""

    def getName(self):
        return self.name

class NewModelStyle(ModelStyle):
    def __init__(self, pdh):
        self.name = 'new'
        self.pdh = pdh

    def getSourcePath(self, basePath, pcssModel):
        """Get model source path; here each modbase sequence has subdirectories for models and alignments in its own sequence directory"""
        return os.path.join(basePath, self.pdh.getThreeLetterOutputDir(pcssModel.getSequenceId()), "models")

class OldModelStyle(ModelStyle):
    def __init__(self, pdh):
        self.name = 'old'
        self.pdh = pdh
                            
class PcssTempDirectory:

    """Class for managing creation of temporary directories and changing back and forth"""

    def __init__(self):
        """Create temporary directory and change into it"""
        self.currentDir = os.getcwd()
        self.tempDir = tempfile.mkdtemp()
        os.chdir(self.tempDir)

    def changeBack(self, removeTemp=True):
        """Change back to original directory; remove temporary by default"""
        os.chdir(self.currentDir)
        if removeTemp:
            shutil.rmtree(self.tempDir)


class PcssDirectoryHandler:

    """Class for all of the directory creation and processing that occurs over the course of a PCSS run"""
    
    def __init__(self, pcssConfig, internalConfig):
        self.pcssConfig = pcssConfig
        self.internalConfig = internalConfig

    def makeTempDirectory(self):
        ptd = PcssTempDirectory()
        return ptd


    def getPcssErrorFile(self):
        
        return self.getFullOutputFile(self.internalConfig["pcss_error_output_file"])

    def getInternalErrorFile(self):
        return self.getFullOutputFile(self.internalConfig["internal_error_output_file"])

    def createOutputDirectory(self):
        """Creates the 'run directory' which is where all output for this program goes. Run directory is created as $run_directory/$run_name
        where each variable specified in the parameter file. Thus if the user wants to run a program twice with different input, simply change
        the run_name parameter and nothing from previous runs will be overwritten
        """
        outputDir = self.pcssConfig["run_directory"]
        runName = self.pcssConfig["run_name"]
        fullOutputDir = os.path.join(outputDir, runName)
        self.fullOutputDir =  os.path.join(outputDir, runName)
        if (not os.path.exists(self.fullOutputDir)):
            os.mkdir(self.fullOutputDir)
            
    def getFullClusterOutputFile(self, fileName):
        return self.getFullOutputFile(fileName)
    
    def getStructureDirectory(self):
        """Return full path of directory where models are stored and copied to"""
        if (not os.path.exists(self.getFullOutputFile("structures"))):
            os.mkdir(self.getFullOutputFile("structures"))
        return self.getFullOutputFile("structures")

    def getFullOutputFile(self, fileName):
        """Return full path fileName where the path is the full run directory for this run"""
        return os.path.join(self.fullOutputDir, fileName)


    def getFullModelFileFromId(self, modelId):
        """Return full path model file name for this modelId"""
        return os.path.join(self.getStructureDirectory(), "%s.pdb" % modelId)

    def getFullModelFile(self, pcssModel):
        """Return the full path model pdb file for this pcssModel"""
        return os.path.join(self.getStructureDirectory(), pcssModel.getPdbFileName())

    def getTwoLetterOutputDir(self, modbaseSequence):
        """Get two letter dirctory whose name is the first two characters of the sequence"""
        twoLetter = modbaseSequence[0:2]
        return twoLetter

    def getThreeLetterOutputDir(self, modbaseSequence):
        """Get directory hierarchy for this sequence where the parent directory is the first three characters of the sequence"""
        threeLetter = modbaseSequence[0:3]
        return os.path.join(threeLetter, modbaseSequence)

    def tryShutil(self, function, *args):
        """Run a python shutil module command"""
        try:
            function(*args)
        except IOError, e:
            raise pcssErrors.PcssShutilError(e, function, args)

    def copyFile(self, sourceDir, sourceFile, destinationDir):
        """Copy source file from sourceDir to destinationDir; don't return until copy is done"""
        
        self.tryShutil(shutil.copy, os.path.join(sourceDir, sourceFile), destinationDir)
        destinationFile = os.path.join(destinationDir, sourceFile)

        #possible error: if file is large, conceivably this could return before copy is done
        self.sleepUntilDone(destinationFile, predicate=self.fileDoesNotExist)

    def moveFile(self, sourceDir, sourceFile, destinationDir):
        """Move source file from sourceDir to destinationDir; don't return until move is done"""
        self.tryShutil(shutil.move, os.path.join(sourceDir, sourceFile), destinationDir)
        destinationFile = os.path.join(destinationDir, sourceFile)
        self.sleepUntilDone(destinationFile, predicate=self.fileDoesNotExist)

    def runSubprocess(self, args, checkStdError=True):
        """Run python subprocess module command; by default, raise exception if anything was written to stderr"""
        process = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
        processOutput = process.communicate()
        if (processOutput[1] != "" and checkStdError):
            raise pcssErrors.PcssGlobalException("Got subprocess error.\nRan method args %s\nGot stderr %s" % (args, processOutput[1])) 
        return processOutput

    def unzipFile(self, sourceFile):
        """Unzip sourceFile; expects .gz suffix"""
        if (not sourceFile.endswith(".gz")):
            #might end up having to change this later to include other filetypes, but this could wreak havoc
            #if result file isn't properly formatted
            raise pcssErrors.PcssGlobalException("Attempted to unzip file %s that does not end with '.gz'")
        resultFile = sourceFile.rstrip(".gz")
        if (not os.path.exists(resultFile)): 
            #should never already exist since we wouldn't be here if it did, but another process could possibly have put it here
            #there is still somewhat of a race condition though as this check could return false and another process could put 
            #the unzipped file in immdediately after, but chances are extremely low
            self.runSubprocess(['gunzip', sourceFile])

            self.sleepUntilDone(resultFile, predicate=self.fileDoesNotExist)

    def fileDoesNotExist(self, fileName):
        """Simple predicate that returns true if fileName doesn't exist in file system"""
        
        return not os.path.exists(fileName)

    def fileExists(self, fileName):
        """Simple predicate that returns true if fileName exists in file system"""
        return os.path.exists(fileName)

    def sleepUntilDone(self, fileName, predicate):
        """Sleep until predicate involving fileName is true; useful for avoiding race conditions in file manipulation"""
        sleepTime = 0
        while (predicate(fileName)):
            print "sleep 1 second"
            time.sleep(1)
            sleepTime += 1
            if (sleepTime > 10):
                raise pcssErrors.PcssGlobalException("Timeout on file %s" % fileName)
     
    def getSvmApplicationSetFile(self):
        return self.getFullOutputFile(self.internalConfig["application_set_file_name"])

    def getSvmTestSetFile(self):
        return self.getFullOutputFile(self.internalConfig["test_set_file_name"])

    def getSvmTrainingSetFile(self):
        return self.getFullOutputFile(self.internalConfig["training_set_file_name"])

    def getSvmNewModelFile(self):
        return self.getFullOutputFile(self.internalConfig["training_new_model_name"])

    def getSvmApplicationOutputFile(self):
        return self.getFullOutputFile(self.internalConfig["application_set_output_file_name"])

    def getFullBenchmarkResultFileName(self):
        return self.getFullOutputFile("%s_%s" % (self.getRunName(), self.internalConfig["benchmark_result_file_suffix"]))

    def getFinalSvmApplicationResultFile(self):
        return self.getFullOutputFile
                      
    def getLeaveOneOutResultFileName(self):
        return self.getFullOutputFile("%s_%s" % (self.getRunName(), self.internalConfig["loo_result_file_suffix"]))

    def getUserModelFileName(self):
        return self.getFullOutputFile("%s_%s" % (self.getRunName(), self.internalConfig["user_model_suffix"]))

    def getSvmTestOutputFile(self):
        return self.getFullOutputFile(self.internalConfig["test_set_output_file_name"])
   
    def getRunName(self):
        return self.pcssConfig["run_name"]
        
    def getSeqBatchDirectory(self):
        directoryName = self.getFullOutputFile(self.internalConfig["seq_batch_directory"])
        if (not os.path.exists(directoryName)):
            os.mkdir(directoryName)
        return directoryName

    def getSeqBatchSubDirectoryName(self, i):
        return os.path.join(self.getSeqBatchDirectory(), str(i))

    def getClusterSeqBatchSubDirectoryName(self, i):
        return getSeqBatchSubDirectoryName(i)

    def getClusterSeqBatchDirectory(self):
        return self.getSeqBatchDirectory()

    def getSeqBatchFastaFileName(self, i):
        subDirectoryName = self.getSeqBatchSubDirectoryName(i)
        fileName = os.path.join(subDirectoryName, self.internalConfig["seq_batch_input_fasta_file_name"])
        return fileName

    def getClusterSeqBatchFastaFileName(self, i):
        return self.getSeqBatchFastaFileName(i)

    def getPcssClusterBaseDirectory(self):
        return self.internalConfig["pcss_directory"]

    def getFullBenchmarkModelFile(self, fileName):
        return os.path.join(self.getPcssClusterBaseDirectory(), "data", "benchmark", fileName)

    def getBenchmarkScoreFile(self):
        return self.pcssConfig["svm_benchmark_file"]

    def getModelFileName(self):
        return self.pcssConfig["svm_model_file"]

    def getClusterNodeConfig(self, i):
        baseConfig = copy.deepcopy(self.pcssConfig)

        baseConfig["pcss_directory"] = self.getPcssClusterBaseDirectory()
        baseConfig["fasta_file"] =  self.getClusterSeqBatchFastaFileName(i)
        
        baseConfig["run_name"] = str(i)
        baseConfig["run_directory"] = self.getFullClusterOutputFile(self.internalConfig["seq_batch_directory"])
        baseConfig["using_web_server"] = False

        return baseConfig



class PcssServerDirectoryHandler(PcssDirectoryHandler):

    def getClusterNodeConfig(self, i):

        baseConfig = copy.deepcopy(self.pcssConfig)
        
        baseConfig["pcss_directory"] = self.getPcssClusterBaseDirectory()
        baseConfig["fasta_file"] =  self.getClusterSeqBatchFastaFileName(i)
        baseConfig["svm_benchmark_file"] = self.getBenchmarkScoreFile()
        baseConfig["svm_model_file"] = self.getModelFileName()
        
        baseConfig["run_name"] = str(i)
        baseConfig["run_directory"] = self.getFullClusterOutputFile(self.internalConfig["seq_batch_directory"])
        baseConfig["using_web_server"] = False

        return baseConfig


    def createOutputDirectory(self):
        self.fullOutputDir = self.pcssConfig["job_directory"]
        #has already been created by job class

    def getClusterSeqBatchDirectory(self):
        seqBatchDir = self.getFullClusterOutputFile(self.internalConfig["seq_batch_directory"])
        return seqBatchDir

    def getClusterSeqBatchFastaFileName(self, i):
        return os.path.join(self.getClusterSeqBatchSubDirectoryName(i), self.internalConfig["seq_batch_input_fasta_file_name"])

    def getClusterSeqBatchSubDirectoryName(self, i):
        return os.path.join(self.getClusterSeqBatchDirectory(), str(i))

    def getFullClusterOutputFile(self, fileName):

        headNodeDir = self.getClusterRunDirectory()
        runName = self.getRunName()
        return os.path.join(headNodeDir, runName, fileName)

    def getClusterRunDirectory(self):
        return self.internalConfig["netapp_server_run_directory"] 

    def getPcssClusterBaseDirectory(self):
        return self.internalConfig["netapp_server_base_directory"]

    def getBenchmarkScoreFile(self):
        bmm = self.setupBenchmarkModelMap()
        frontendName = self.pcssConfig["svm_application_model"]    
        return bmm.getBenchmarkScoreFile(frontendName)

    def getModelFileName(self):
        bmm = self.setupBenchmarkModelMap()
        frontendName = self.pcssConfig["svm_application_model"]    
        return bmm.getModelFileName(frontendName)

    def setupBenchmarkModelMap(self):
        if "using_custom_model" in self.pcssConfig:
            if (self.pcssConfig["using_custom_model"] == "False"):
                bmm = BenchmarkModelMap(self)
                return bmm

class BenchmarkModelMap:
    def __init__(self, pdh):
        self.pdh = pdh
        fileReader = PcssFileReader(self.pdh.internalConfig["benchmark_model_map_file_name"])
        lines = fileReader.getLines()
        self.ModelMapTuple = myCollections.namedtuple('modelMap', ['frontendName', 'benchmarkScoreName', 'modelFileName'])
        self.modelMap = {}
        for line in lines:
            [frontendName, benchmarkScoreName, modelFileName] = line.split('\t')
            mmt = self.ModelMapTuple(frontendName, benchmarkScoreName, modelFileName)
            self.modelMap[frontendName] = mmt

    def getBenchmarkScoreFile(self, frontendName):
        benchmarkScoreFile = self.modelMap[frontendName].benchmarkScoreName
        return self.pdh.getFullBenchmarkModelFile(benchmarkScoreFile)

    def getModelFileName(self, frontendName):
        modelFileName = self.modelMap[frontendName].modelFileName
        return self.pdh.getFullBenchmarkModelFile(modelFileName)


class PcssFileReader:

    """Essentially python file object that provides commonly used processing functionality"""

    def __init__(self, fileName, skipBlanks=True, skipHash=True):
        """Read input file and save all lines after processing

        Similar to python File object, except we commonly will want to skip blank lines in the file
        (programs often break when getting an unexpected blank line) and also skip lines beginning with
        '#' as these are often comments as oposed to data.

        @param fileName full path file name to be processed
        @param skipBlanks set to true to ignore blank lines
        @param skipHash set to true to ignore lines beginning with '#'
        """
        fileHandle = open(fileName, 'r')
        lines = fileHandle.readlines()
        finalLines = []
        blankRe = re.compile('^\s*$')
        hashRe = re.compile('^\#')

        for line in lines:
            line = line.rstrip('\r\n')
            if (self.skipLine(line, blankRe, skipBlanks)):
                continue
            if (self.skipLine(line, hashRe, skipHash)):
                continue
            
            finalLines.append(line)
        self.finalLines = finalLines
        self.fileName = fileName
        fileHandle.close()
    def skipLine(self, line, regex, skip):
        """Return True if this line should be skipped in processing"""
        testLine = regex.search(line)
        if (testLine and skip):
            return True
        return False
        

    def getLines(self):
        return self.finalLines

class FastaGrabber:

    def readSourceFastaFile(self):
        fh = open(self.sourceFastaFile, 'r')
        self.fastaIterator = SeqIO.FastaIO.FastaIterator(fh)

    def getSeqList(self, sequence, size):
        if (size == 0):
            return [sequence]
        return split_len(sequence, size)

class TwoLetterFastaGrabber(FastaGrabber):
    
    def __init__(self, sourceFastaFile, twoLetterCode, outputFileName, size=60):
        
        self.sourceFastaFile = sourceFastaFile
        self.readSourceFastaFile()
        outputFh = open(outputFileName, 'w')
        seqCount = 0
        for seqRecord in self.fastaIterator:
            if (seqRecord.id.startswith(twoLetterCode)):
                seqCount += 1
                outputFh.write(">%s\n" % seqRecord.id)
                seqList = self.getSeqList(str(seqRecord.seq), size)
                outputFh.write("%s\n" % '\n'.join(seqList))
        outputFh.close()
        print "wrote %s sequences with code %s to output file %s" % (seqCount, twoLetterCode, outputFileName)
