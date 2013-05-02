import sys
import os
from Bio import SeqIO
from Bio import PDB
import itertools
import pcssTools
import StringIO
from operator import itemgetter
import logging
import subprocess
import pcssErrors
import tempfile
import shutil
log = logging.getLogger("pcssModels")



class PcssModelTableColumns:
    
    """Class handling order of column names that are read in from the model table; helps provide controlled vocabulary for subsequent access"""

    def __init__(self, pcssConfig):
        reader = pcssTools.PcssFileReader(pcssConfig['model_table_column_file'])
        lines = reader.getLines()
        columnOrder = 0
        self._columnDict = {}
        for line in lines:
            self._columnDict[columnOrder] = line
            columnOrder += 1

    def getColumnNameIterator(self):
        """Get column names sorted by how they appeared in the input column order file"""
        columnNames = []
        for columnNumber, columnName in (sorted(self._columnDict.iteritems(), key=itemgetter(0))):
            columnNames.append(columnName)
        return columnNames

    def getColumnName(self, columnNumber):
        return self._columnDict[columnNumber]

    def getColumnCount(self):
        return len(self._columnDict.keys())

class PcssModelTable:

    """Class representing a model table file that was generated from an SQL query; creates models and groups them by their proteins"""

    def __init__(self, pcssRunner, modelTableColumns):
        """Read model table file and make one model for each line"""
        reader = pcssTools.PcssFileReader(pcssRunner.pcssConfig['model_table_file'])
        lines = reader.getLines()
        self.pcssRunner = pcssRunner
        self._sequenceDict = {}
        for line in lines:
            pcssModel = PcssModel(pcssRunner)
            pcssModel.initFromModelTableLine(line, modelTableColumns)
            
            self.addModel(pcssModel)

    def addModel(self, pcssModel):
        sequenceId = pcssModel.getAttributeValue("seq_id")
        modelSequence = self.getPcssModelSequence(sequenceId)
        modelSequence.addModel(pcssModel)

    def getPcssModelSequence(self, sequenceId):
        """Factory method for getting a PcssModelSequence for the input sequenceId"""
        modelSequence = None

        if (sequenceId in self._sequenceDict):
            return self._sequenceDict[sequenceId]
        else:
            modelSequence = PcssModelSequence(sequenceId)
            self._sequenceDict[sequenceId] = modelSequence
        return modelSequence
            
    def getSequences(self):
        return self._sequenceDict.values()

class PcssModelSequence:

    """Simple class for grouping models read in by model table; models from here are later transferred to PcssProteins"""

    def __init__(self, modbaseSeqId):
        self._models = []
        self.modbaseSeqId = modbaseSeqId

    def addModel(self, pcssModel):
        self._models.append(pcssModel)

    def getModels(self):
        return self._models

    def getModel(self, modelId):
        for m in self._models:
            if (m.getId() == modelId):
                return m

    def getOutput(self):
        result = "\nModbase Seq: %s\n" % self.modbaseSeqId
        modelOutputList = []
        for model in self._models:
            modelOutputList.append(model.getOutput())
            
        result += "\n".join(modelOutputList)
        return result

class PcssModel:

    """Class for a homology model for a protein. Stores all attributes for the model and handles DSSP processing"""

    def __init__(self, pcssRunner):
        self._attributes = {}
        self.pcssRunner = pcssRunner
        self.bioModel = None
        self.dssp = None
        
    def loadDsspResults(self):
        """Get the PDB file for this model and use it as input for DSSP"""
        self.loadBioModelPdb()
        self.runDssp()

    def runDssp(self):
        """Run DSSP executable for this model"""
        if (self.dssp is None):
            dssp = PDB.DSSP(self.bioModel, self.pcssRunner.pdh.getFullModelFile(self), 
                            self.pcssRunner.pcssConfig["dssp_executable"])
            #Hard to get exact reason why DSSP didn't work since it's BioPython, but will set 
            #as feature exception rather than global exception
            if (dssp is None or len(dssp.keys()) < 1):
                raise pcssErrors.DsspException("Did not load DSSP for model %s. This likely indicates a problem with the "
                                               "Biopython DSSP module;\ntry running DSSP from the command line to isolate "
                                               "the issue" % self.getId())
            self.dssp = dssp

    def getRelativeSolventAcc(self, residueIndex):
        """Return the fraction of the residue that is accessible to solvent according to DSSP"""
        dsspList = self.getDsspTuple(residueIndex)
        relativeSolventAcc = dsspList[3]
        return relativeSolventAcc

    def getSecondaryStructure(self, residueIndex):
        """Return the secondary structure call for this residue according to DSSP"""
        dsspList = self.getDsspTuple(residueIndex)
        secondaryStructure = dsspList[1]
        return secondaryStructure

    def getDsspResidueCode(self, residueIndex):
        """Get the one-letter residue code for this residue according to DSSP (mostly for testing / validation"""
        dsspList = self.getDsspTuple(residueIndex)
        residueCode = pcssTools.getOneLetterFromBioResidue(dsspList[0].get_resname())

        return residueCode

    def getDsspTuple(self, residueIndex):
        """Return DSSP information for this residue.

        DSSP is sort of messy; use this as a convenience method for accessing its results. Results are returned as a 
        5-element list: [residue object, secondary structure, solvent accessibility, relative solvent accessibility, ?]"""
        dsspKey = (' ', self.bioChain._translate_id(residueIndex + 1))
        dsspTuple =  self.dssp[dsspKey]
        return dsspTuple

    def loadBioModelPdb(self):
        """Create a BioPython model and chain object for this model"""
        if (self.bioModel is None):
            modelFileName = self.pcssRunner.modelHandler.getLocalModelFileName(self)
            parser =  self.pcssRunner.getBioParser()
            self.bioModel = parser.get_structure(self.getId(), modelFileName)[0]
            self.bioChain =  self.bioModel[" "]
                
    def getPdbFileName(self):
        return "%s.pdb" % self.getId()

    def getRunName(self):
        return self.getAttributeValue("run")

    def getSequenceId(self):
        return self.getAttributeValue("seq_id")

    def getAttributeValue(self, attributeName):
        """"Return string for the given attribute as it appears in the model table"""
        if (not attributeName in self._attributes):
            raise pcssErrors.PcssGlobalException("Error: model does not have attribute %s. Make sure this attribute "
                                                 "is specified in the model column info file" % attributeName)
        return self._attributes[attributeName]


    def getAttributeNames(self):
        return self._attributes.keys()

    def setAttribute(self, attributeName, attributeValue):
        self._attributes[attributeName] = attributeValue

    def initFromModelTableLine(self, line, modelTableColumns):
        """Initialize model from file.

        Use ModelTableColumns to get the names of each attribute in the line (which is read from the model table file)
        and save these attributes internally"""
        cols = line.split('\t')
        i = 0

        if (modelTableColumns.getColumnCount() != len(cols)):
            raise pcssErrors.PcssGlobalException("Model table column order file contains a "
                                                 "different number of columns (%s) than model table (%s)\nline: %s" % 
                                                 (modelTableColumns.getColumnCount(), len(cols), line))
        for col in cols:
            columnName = modelTableColumns.getColumnName(i)
            self.setAttribute(columnName, col)
            i += 1

    def containsPeptide(self, peptide):
        """Return true if the given peptide is fully contained within this model"""
        return (peptide.startPosition > self.getModelStart() and
                peptide.endPosition < self.getModelEnd())
    
    def getModelStart(self):
        return int(self.getAttributeValue("target_beg"))
    
    def getModelEnd(self):
        return int(self.getAttributeValue("target_end"))

    def getRange(self):
        return [self.getAttributeValue("target_beg"), self.getAttributeValue("target_end")]

    def getId(self):
        return self.getAttributeValue("model_id")

    def getOutput(self):
        outputList = []
        for attributeName, attributeValue in self._attributes.iteritems():
            outputList.append("%s = %s" % (attributeName, attributeValue))

        return "; ".join(outputList)

    def calculateCoverage(self, proteinSeqLength):
        """Calculate the fraction of thep protein that is covered by this model"""
        self.setAttribute("coverage", float(self.getLength()) / float(proteinSeqLength))
    
    def setModelUrl(self, modelUrl):
        self.setAttribute("model_url", modelUrl)

    def getLength(self):
        targetEnd = int(self.getAttributeValue("target_end"))
        targetBegin = int(self.getAttributeValue("target_beg"))
        if (targetEnd <= targetBegin):
            raise pcssErrors.PcssGlobalException("Error in model table: model %s has target_end position %s "
                                                 "before target_start position %s" % (self.getId(), targetEnd, targetBegin))
        return targetEnd - targetBegin 

    def isEqual(self, otherModel):
        if (otherModel is None):
            return False
        for attName in otherModel.getAttributeNames():
            if (attName not in self._attributes):
                return False
            if (str(self.getAttributeValue(attName)) != str(otherModel.getAttributeValue(attName))):
                return False
        return True
