import sys
from validate import ValidateError
import os
import re
from Bio import SeqIO
import pcssIO
import pcssCluster
import tempfile
import configobj
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
#from pympler import summary
#from pympler import muppy
#from pympler import tracker
import traceback
from Bio import PDB
log = logging.getLogger("pcssTools")

#tr = tracker.SummaryTracker()
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
        print protein.getAttributeOutputString("protein_errors")
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
        self.internalConfig = configobj.ConfigObj(pcssConfig["internal_config_file"], configspec=pcssConfig["internal_config_spec_file"])
        self.validateConfig(self.internalConfig)

        self.pdh = PcssDirectoryHandler(pcssConfig, self.internalConfig)
        self.pdh.createOutputDirectory()
        #self.modelHandler = PcssModelHandler(pcssConfig, self.pdh)
        self.parser = PDB.PDBParser(QUIET=True)
        self.pfa = pcssIO.PcssFileAttributes(pcssConfig)

        self.initSubclass()

    def initSubclass(self):
        return

    def execute(self):
        self.checkForErrors()
        try:
            self.executePipeline()
        except pcssErrors.PcssGlobalException as pge:
            print pge.msg
            tb = traceback.format_exc()
            print "WRITING EXCEPTION " + pge.msg + "\n" + tb
            self.writePcssErrorFile(pge.msg + "\n" + tb)

        except pcssErrors.ErrorExistsException:
            return

        except pcssErrors.InternalException as e:
            self.writeInternalErrorFile(e.msg)
        
        except Exception as e:
            print "WRITING INTERNAL EXCPETION"
            print e
            self.writeInternalErrorFile(e)


    def validateConfig(self, config):
        validator = Validator({'file': fileExists})
        results = config.validate(validator, preserve_errors=True)
        if (results != True):
            self.handleConfigError(results)
            
    def checkForErrors(self):
        if (os.path.exists(self.pdh.getPcssErrorFile()) or 
            os.path.exists(self.pdh.getInternalErrorFile())):
            raise pcssErrors.ErrorExistsException("Previous error exists")


    def validatePeptideCodeStatus(self, status, peptideCode):
        try:
            return self.validatePeptideStatus(status)
        except pcssErrors.PcssGlobalException as e:
            raise pcssErrors.PcssGlobalException("%s (peptide code: %s)" % (e.msg, peptideCode))

    def validatePeptideStatus(self, status):
                                                
        status = status.lower()
        if (not (status == self.getPositiveKeyword() or status == self.getNegativeKeyword() or status == self.getApplicationKeyword())):
            raise pcssErrors.PcssGlobalException("Peptide status %s not valid status (needs to be %s, %s, %s" % (status, self.getPositiveKeyword(), 
                                                                                                                 self.getNegativeKeyword, self.getApplicationKeyword))
        return status

    def validatePeptideTrainingStatus(self, status):
                                                
        status = status.lower()
        if (not (status == self.getPositiveKeyword() or status == self.getNegativeKeyword())):
            raise pcssErrors.PcssGlobalException("Peptide status %s not valid status (needs to be %s, %s)" % (status, self.getPositiveKeyword(), 
                                                                                                              self.getNegativeKeyword()))
        return status

    def getBioParser(self):
        return self.parser

    def getSvmFeatureOrder(self):
        featureOrderList = self.internalConfig["feature_order"]
        return featureOrderList

    def writePcssErrorFile(self, message):
        errorFile = self.pdh.getPcssErrorFile()
        errorFh = open(errorFile, 'w')
        errorFh.write(self.internalConfig["keyword_pcss_error"] + "\n")
        errorFh.write(message + "\n")
        tb = traceback.format_exc()
        errorFh.write(tb + "\n")

    def writeInternalErrorFile(self, e):
        errorFile = self.pdh.getInternalErrorFile()
        errorFh = open(errorFile, 'w')
        errorFh.write(self.internalConfig["keyword_internal_error"] + "\n")
        errorFh.write(str(e) + "\n")
        tb = traceback.format_exc()
        errorFh.write(tb + "\n")

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

    def getSeqBatchDirectory(self):
        directoryName = self.pdh.getFullOutputFile(self.internalConfig["seq_batch_directory"])
        if (not os.path.exists(directoryName)):
            os.mkdir(directoryName)
        return directoryName
        
    def getJobDirectory(self):
        return self.pdh.getJobDirectory()

    def setJobDirectory(self, directoryName):
        self.pdh.setJobDirectory(directoryName)

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
        msg = ""
        print "An error occurred when validating a configuration file. Please check the log file for more details."
        for (section_list, key, _) in flatten_errors(self.pcssConfig, results):
            if key is not None:
                msg +=  'The "%s" key in the section "%s" failed validation\n' % (key, ', '.join(section_list))
            else:
                msg += 'The following section was missing:%s ' % ', '.join(section_list)
        print msg
        raise pcssErrors.PcssGlobalException(msg)
                
    def addPeptideFeatures(self):
        
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
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

class DisopredStandaloneRunner(PcssRunner):
    
    def executePipeline(self):

        self.readProteins()

        self.runDisopred()

    def runDisopred(self):

        disopredFileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.pdh)
        disopredReader = pcssFeatureHandlers.DisopredReader(disopredFileHandler)
        disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(disopredFileHandler)

        for protein in self.proteins:
            log.debug("begin processing disopred for protein %s" % protein.modbaseSequenceId)
            protein.processDisopred(disopredReader, disopredRunner)
            log.debug("done processing disopred for protein %s" % protein.modbaseSequenceId)

class ModelRunner(PcssRunner):
    def initSubclass(self):
        self.modelHandler = PcssModelHandler(self.pcssConfig, self.pdh)

class SvmApplicationClusterRunner(PcssRunner):
    def executePipeline(self):

        seqDivider = pcssCluster.SeqDivider(self)
        
        seqDivider.divideSeqsFromFasta()
        
        seqDivider.makeFullSvmApplicationSgeScript()

class TrainingBenchmarkClusterRunner(PcssRunner):

    def executePipeline(self):
        tcb =  pcssCluster.ClusterBenchmarker(self)
            
        tcb.prepareTrainingBenchmarkRun()
            
        tcb.makeFullTrainingBenchmarkScript()
        
class TrainingAnnotationClusterRunner(PcssRunner):
    def executePipeline(self):
        seqDivider = pcssCluster.SeqDivider(self)
            
        seqDivider.divideSeqsFromFasta()
            
        seqDivider.makeFullTrainingAnnotationSgeScript()
            
class FinalizeApplicationClusterRunner(PcssRunner):
    def executePipeline(self):
        seqDivider = pcssCluster.SeqDivider(self)

        seqDivider.mergeSvmApplicationResults()
        
class SvmApplicationRunner(ModelRunner):
    def runSvm(self):
        self.appSvm = pcssSvm.ApplicationSvm(self)
        self.appSvm.setProteins(self.proteins)
        self.appSvm.writeClassificationFile()
        self.appSvm.classifySvm()
        self.appSvm.readResultFile()
        self.appSvm.addScoresToPeptides()

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

class CompleteSvmRunner(PcssRunner):
    def executePipeline(self):
        self.readAnnotationFile()
        
        self.createSvmModel()

    def createSvmModel(self):
        generator = pcssSvm.CompleteSvmGenerator(self)

        generator.createSvmModel(getAllPeptides(self.reader.getProteins(), False))

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

    
class TrainingBenchmarkRunner(PcssRunner):
    def executePipeline(self):

        self.readAnnotationFile()
        
        self.benchmark()
              
    def benchmark(self):
        
        benchmarker = pcssSvm.SvmBenchmarker(self)
 
        for i in range(self.pcssConfig["training_iterations"]):
            benchmarker.createTrainingAndTestSets(getAllPeptides(self.proteins, False))  #parition pepitdes, write training and test set files

            benchmarker.trainAndApplyModel() #call svm command line

            benchmarker.readBenchmarkResults() #read result file

        benchmarker.processAllResults()

class PcssModelHandler:

    """Class for managing model PDB files, retrieving them from file servers as necessary"""

    def __init__(self, pcssConfig, pdh):
        self.pcssConfig = pcssConfig
        self.pdh = pdh
        self.loadRunInfo()


    def loadRunInfo(self):
        """Load PcssModelRunInfo object which provide data specific to differen modpipe runs"""
        
        self._runInfo = {}
        runInfoFile  = self.pcssConfig["model_run_info"]
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
        if (not self.pcssConfig["using_web_server"]):  #consider changing to 'head_node'
            outputDir = self.pcssConfig["run_directory"]
            runName = self.pcssConfig["run_name"]
            fullOutputDir = os.path.join(outputDir, runName)
            self.fullOutputDir =  os.path.join(outputDir, runName)
            if (not os.path.exists(self.fullOutputDir)):
                os.mkdir(self.fullOutputDir)
        else:
            self.fullOutputDir = self.pcssConfig["job_directory"]

            
    def getJobDirectory(self):
        return self.jobDirectory

    def setJobDirectory(self, directoryName):
        self.jobDirectory = directoryName
        if (not os.path.exists(directoryName)):
            os.mkdir(directoryName)

    def getStructureDirectory(self):
        """Return full path of directory where models are stored and copied to"""
        return self.pcssConfig["model_directory"]

    def getFullOutputFile(self, fileName):
        """Return full path fileName where the path is the full run directory for this run"""
        if (not self.pcssConfig["using_web_server"]):  #consider changing to 'head_node'
            return os.path.join(self.fullOutputDir, fileName)
        else:
            return os.path.join(self.jobDirectory, fileName)

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
        except IOError as e:
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
        return self.getFullOutputFile("%s_%s" % (self.pcssConfig["run_name"], self.internalConfig["benchmark_result_file_suffix"]))

    def getFinalSvmApplicationResultFile(self):
        return self.getFullOutputFile
                      
    def getLeaveOneOutResultFileName(self):
        return self.getFullOutputFile("%s_%s" % (self.pcssConfig["run_name"], self.internalConfig["loo_result_file_suffix"]))

    def getUserModelFileName(self):
        return self.getFullOutputFile("%s_%s" % (self.pcssConfig["run_name"], self.internalConfig["user_model_suffix"]))

    def getSvmTestOutputFile(self):
        return self.getFullOutputFile(self.internalConfig["test_set_output_file_name"])
   
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
