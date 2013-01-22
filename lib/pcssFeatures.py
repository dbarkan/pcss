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
import pcssSvm
import tempfile
import shutil
log = logging.getLogger("pcssFeatures")

class SequenceFeatureCallSet:
    """Simple class for storing a series of sequence features (psipred / disopred) on individual residues"""
    def __init__(self):
        self._sequenceFeatureCalls= {}

    def addSequenceFeatureCall(self, sequenceFeatureCall):
        """Add a call for an individual residue to my set"""
        self._sequenceFeatureCalls[sequenceFeatureCall.residueNumber] = sequenceFeatureCall

    def getSequenceFeatureCall(self, baseOneResidueNumber):
        if (not self._sequenceFeatureCalls.has_key(baseOneResidueNumber)):
            raise self.getFeatureNotFoundException(baseOneResidueNumber)
        return self._sequenceFeatureCalls[baseOneResidueNumber]
    
    def checkSequenceMatch(self, proteinSequence):
        """Check over all of my calls whether each has the same amino acid as expected in the input protein sequence"""
        for sequenceFeatureCall in self._sequenceFeatureCalls.values():
            sequenceFeatureCall.checkSequenceMatch(proteinSequence)

    def makeFullCallString(self):
       """Return all calls as a string of one letter call codes"""
       callString = ""
       for residueNumber, sequenceFeatureCall in (sorted(self._sequenceFeatureCalls.iteritems(), key=itemgetter(0))):
           callString += sequenceFeatureCall.call
       return callString

class DisopredSequenceFeatureCallSet(SequenceFeatureCallSet):
    def getFeatureNotFoundException(self, startResidue):
        return pcssErrors.DisopredPeptideNotFoundException("Disopred result file did not contain peptide start residue %s" % startResidue)

class PsipredSequenceFeatureCallSet(SequenceFeatureCallSet):
    def getFeatureNotFoundException(self, startResidue):
        return pcssErrors.PsipredPeptideNotFoundException("Psipred result file did not contain peptide start residue %s" % startResidue)


class SequenceFeatureCall:
    def checkSequenceMatch(self, proteinSequence):
        """Perform sanity check to make sure protein sequence residue matches mine at the same position"""
        proteinSequenceAA = proteinSequence[self.residueNumber - 1]
        if (self.residueOneLetter != proteinSequenceAA):
            raise self.getMismatchException("%s mismatch between result file (%s) and protein sequence (%s) at position %s" % 
                                            (self.name, self.residueOneLetter, proteinSequenceAA, self.residueNumber))                                           

class DisorderResidueCall(SequenceFeatureCall):

    """Simple class for storing a Disopred call for one residue"""
    
    def __init__(self, disopredLine):
        """@param disopredLine: one line from the disopred result file"""

        if (len(disopredLine.split()) != 5):
            raise pcssErrors.DisopredBadLineException("Did not get 5 columns out of line %s " % disopredLine)
        [residueNumber, self.residueOneLetter, disorderCall, firstScore, secondScore] = disopredLine.split()
        self.name = "Disopred"
        self.residueNumber = int(residueNumber)
        self.score = float(firstScore)
        if (disorderCall == "*"):
            self.call = "D"
        elif (disorderCall == "."):
            self.call = "O"
        else:
            raise pcssErrors.DisopredBadCallException("Disopred got unexpected call of %s from line %s (expect * or .)" % (disorderCall, disopredLine))

    def getMismatchException(self, msg):
        return pcssErrors.DisopredMismatchException(msg)

class PsipredResidueCall(SequenceFeatureCall):

    """Simple class for storing a Psipred call for one residue"""

    def __init__(self, psipredLine):
        """@param psipredLine: one line from the psipred result file"""
        if (len(psipredLine.split()) != 6):
            raise pcssErrors.PsipredBadLineException("Did not get 6 columns out of line %s " % psipredLine)

        [residueNumber, self.residueOneLetter, psipredCall, firstScore, secondScore, thirdScore] = psipredLine.split()
        self.name = "Psipred"
        self.residueNumber = int(residueNumber)
        self.score = float(firstScore)
        if (psipredCall == "C"):
            self.call = "L"
        elif (psipredCall == "H"):
            self.call = "A"
        elif (psipredCall == "E"):
            self.call = "B"
        else:
            raise pcssErrors.PsipredBadCallException("Psipred got unexpected call of %s from line %s (expect H,C,E)" % (psipredCall, psipredLine))

    def getMismatchException(self, msg):
        return pcssErrors.PsipredMismatchException(msg)


class PcssFeature:

    """Represents one feature which can just be annotation or can be used for SVM Model input"""

    def getOutputString(self):
        return "%s: %s" % (self.name, self.getValueString())

    def initFromFileValue(self, fileValue):
        print "%s: inititng from %s" % (self.name, fileValue)

    def getEmptyFeatureOffset(self, peptideLength):
        return self.getFeatureLength() * peptideLength

    def convertStringListToFloat(self, stringList):
        floatList = []
        for nextString in stringList:
            floatList.append(float(nextString))
        return floatList

    def isInitialized(self):
        return True

    def getFeatureLength(self):
        return 1

class DisorderStringFeature(PcssFeature):
    def __init__(self, disorderStringList=None):
        self.disorderStringList = disorderStringList
        self.name = "disopred_string_feature"
    
    def getValueString(self):
        if (self.disorderStringList is None):
            return ""
        return "".join(self.disorderStringList)

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.disorderStringList = list(fileValue)
    
class DisorderScoreFeature(PcssFeature):
    def __init__(self, disorderScoreList=None):
        self.name = "disopred_score_feature"
        self.disorderScoreList = disorderScoreList

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.disorderScoreList = self.convertStringListToFloat(fileValue.split(', '))

    def getValueString(self):
        if (self.disorderScoreList is None):
            return ""

        return ", ".join(str(x) for x in self.disorderScoreList)

    def makeSvmFeature(self, svmHandler):
        svmFeatureList = []
        for score in self.disorderScoreList:
            svmFeatureList.append("%s:%s" % (svmHandler.getFeatureNumber(), score))
            svmHandler.processFeature(1)

        return ' '.join(svmFeatureList)

    def isInitialized(self):
        return self.disorderScoreList is not None

class PsipredStringFeature(PcssFeature):
    def __init__(self, psipredStringList=None):
        self.psipredStringList = psipredStringList
        self.name = "psipred_string_feature"
    
    def getValueString(self):
        if (self.psipredStringList is None):
            return ""
        return "".join(self.psipredStringList)
    
    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.psipredStringList = list(fileValue)



class PeptideSequenceFeature(PcssFeature):
    def __init__(self, sequence=None):
        self.sequence = sequence
        self.name = "peptide_sequence"
        self.seqList = []
        if (sequence is not None):
            self.populateSeqList()
        
        self.makeSvmMap()

    def getFeatureLength(self):
        return 20

    def populateSeqList(self):
        for i in range(len(self.sequence)):
            self.seqList.append(self.sequence[i])

    def initFromFileValue(self, fileValue):
        self.sequence = fileValue
        self.populateSeqList()

    def makeSvmMap(self):
        self.residueOrder = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def getSingleResidueFeatureList(self, residueCode, featureNumber, seqList):
        featureList = []
        foundResidue = False
        for (i, nextResidueCode) in enumerate(self.residueOrder):
            if (nextResidueCode == residueCode):
                foundResidue = True
                featureList.append("%s:%s" % (featureNumber + i, 1))
            else:
                featureList.append("%s:%s" % (featureNumber + i, 0))
        if (not foundResidue):
            raise pcssErrors.PcssGlobalException("Residue %s in sequence %s is not one of the 20 standard amino acids" 
                                                 % (residueCode, "".join(seqList)))
        return " ".join(featureList)

    def getValueString(self):
        if (self.seqList is None):
            return ""
        return "".join(self.seqList)

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.seqList = list(fileValue)

    def makeSvmFeature(self, svmHandler):
        featureList = []
        for residue in self.seqList:
            featureList.append(self.getSingleResidueFeatureList(residue, svmHandler.getFeatureNumber(), self.seqList))
            svmHandler.processFeature(20)
        return " ".join(featureList)

class StringAttribute(PcssFeature):
    def __init__(self, name=None, value=None):
        self.setValues(name, value)

    def setValues(self, name, value):
        self.name = name
        self.value = value
    
    def getValueString(self):
        
        return self.value

    def initFromFileValue(self, name, value):
        self.setValues(name, value)

    def isInitialized(self):
        return self.value is not None

    def makeSvmFeature(self):
        if (self.value.startswith(pcssTools.getPeptideErrorCodePrefix())):
            return 

class PsipredScoreFeature(PcssFeature):
    def __init__(self, psipredScoreList=None):
        
        self.psipredScoreList = psipredScoreList
        self.name = "psipred_score_feature"

    def getValueString(self):
        if (self.psipredScoreList is None):
            return ""
        return ", ".join(str(x) for x in self.psipredScoreList)

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.psipredScoreList = self.convertStringListToFloat(fileValue.split(", "))

    def makeSvmFeature(self, svmHandler):
        svmFeatureList = []
        for score in self.psipredScoreList:
            svmFeatureList.append("%s:%s" % (svmHandler.getFeatureNumber(), score))
            svmHandler.processFeature(1)

        return ' '.join(svmFeatureList)
                                 
    def isInitialized(self):
        return self.psipredScoreList is not None

class DsspStructureFeature(PcssFeature):
    def __init__(self, inputDsspStructureList=None):
        self.createStructureMap()
        
        if (inputDsspStructureList is not None):
            self.dsspStructureList = []
            
            for dsspCall in inputDsspStructureList:
                
                mappedCall = self.getMappedCall(dsspCall)
                self.dsspStructureList.append(mappedCall)
        else:
            self.dsspStructureList = None
        self.name = "dssp_structure"

    def createStructureMap(self):
        """DSSP has multiple representations for each structure type (e.g. 3-10 helix vs normal helix); collapse them to one type"""
        self._structureMap = {}
        self._callMap = {}
        self._structureMap["H"] = "A"
        self._structureMap["G"] = "A"
        self._structureMap["I"] = "A"

        self._structureMap["B"] = "B"
        self._structureMap["E"] = "B"
        
        self._structureMap["T"] = "L"
        self._structureMap["S"] = "L"
        self._structureMap["-"] = "L"

        self._callMap["B"] = 3
        self._callMap["A"] = 2
        self._callMap["L"] = 1

    def getMappedCall(self, dsspCall):
        if (dsspCall not in self._structureMap):
            
            raise pcssErrors.DsspException("Error: could not map call from dssp call %s" % dsspCall)
        return self._structureMap[dsspCall]
                 
    def getValueString(self):
        if (self.dsspStructureList is None):
            return ""
        return "".join(self.dsspStructureList)

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.dsspStructureList = list(fileValue)

    def makeSvmFeature(self, svmHandler):
        svmFeatureList = []
        for call in self.dsspStructureList:
            value = self.getValueForCall(call)
            svmFeatureList.append("%s:%s" % (svmHandler.getFeatureNumber(), value))
            svmHandler.processFeature(1)
            
        return ' '.join(svmFeatureList)

    def getValueForCall(self, call):
        if (call not in self._callMap):
            raise pcssErrors.GlobalException("Error: did not get internal mapped value for DSSP structure call %s" % call)
        return self._callMap[call]

    def isInitialized(self):
        return self.dsspStructureList is not None


class PeptideErrorFeature(PcssFeature):

    """An error is its own feature; track errors for peptides using this class"""

    def __init__(self):
        self.errorCodeList = []
        self.name = "peptide_errors"
    
    def getValueString(self):
        
        if (len(self.errorCodeList) == 0):
            return "none"
        else:
            return '; '.join(self.errorCodeList)
    
    def addError(self, errorCode):
        self.errorCodeList.append(errorCode)

    def initFromFileValue(self, fileValue):
        if (fileValue != "none"):
            self.errorCodeList = fileValue.split('; ')

class DsspAccFeature(PcssFeature):
    def __init__(self, dsspAccList=None):
        self.dsspAccList = dsspAccList
        self.name = "dssp_accessibility"

    def initFromFileValue(self, fileValue):
        if (fileValue != ""):
            self.dsspAccList = self.convertStringListToFloat(fileValue.split(', '))

    def getValueString(self):
        if (self.dsspAccList is None):
            return ""
        return ", ".join(str(round(x, 3)) for x in self.dsspAccList)

    def makeSvmFeature(self, svmHandler):

        svmFeatureList = []
        for score in self.dsspAccList:
            svmFeatureList.append("%s:%s" % (svmHandler.getFeatureNumber(), score))
            svmHandler.processFeature(1)
        return ' '.join(svmFeatureList)

    def isInitialized(self):
        return self.dsspAccList is not None

