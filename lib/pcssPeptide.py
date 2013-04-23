import sys
import os
from Bio import SeqIO
import itertools
import pcssTools
import logging
import pcssSvm
import pcssErrors
import pcssFeatures
import pcssFeatureHandlers
log = logging.getLogger("pcssPeptide")

class PcssProtein(object):

    """Class for one protein sequence, containing a number of peptides on which features are assessed"""
    
    def __init__(self, modbaseSequenceId, pcssRunner):
        self.modbaseSequenceId = modbaseSequenceId
        self.peptides = {}
        self._pcssModels = None
        self.uniprotId = "<uninitialized>"
        self.pcssRunner = pcssRunner
        self.proteinSequence = None
        self.proteinAttributes = {}
        self.setStringAttribute("seq_id", modbaseSequenceId)
        self.setStringAttribute("uniprot_id", self.uniprotId)
        self.setStringAttribute("protein_errors", "none")

    def getAttributeOutputString(self, attributeName):
        """Return attribute value as a string"""
        if (not self.hasAttribute(attributeName)):
            return None
        else:
            return self.proteinAttributes[attributeName].getValueString()

    def hasAttribute(self, attributeName):
        return attributeName in self.proteinAttributes

    def hasModels(self):
        return len(self._pcssModels.values()) > 0

    def setStringAttribute(self, attributeName, attributeValue):
        self.pcssRunner.pfa.validateAttribute(attributeName)
        self.proteinAttributes[attributeName] = pcssFeatures.StringAttribute(attributeName, attributeValue)

    def setProteinSequence(self, sequence):
        self.proteinSequence = sequence
        self.setStringAttribute("protein_length", self.getSequenceLength())

    def getSequenceLength(self):
        return len(self.proteinSequence)

    def setPeptides(self, peptideList):
        if (len(peptideList) == 0):
            self.setStringAttribute("protein_errors", "no_peptides_parsed")
        else:
            for nextPeptide in peptideList:
                self.setPeptide(nextPeptide)

    def setPeptide(self, peptide):
        self.peptides[peptide.startPosition] = peptide

    def setUniprotId(self, uniprotId):
        self.uniprotId = uniprotId
        self.setStringAttribute("uniprot_id", uniprotId)

    def getOutputString(self):
        outputString = "Modbase Seq Id: %s " % self.modbaseSequenceId
        outputString += "Uniprot Id: %s " % self.uniprotId
        outputString += "Peptides:\n%s\n" % self.getPeptideListOutputString()
        return outputString

    def hasErrors(self):
        if (self.getAttributeOutputString("protein_errors") == "none"):
            return False
        return True
    
    def hasPeptide(self, startPosition):
        return startPosition in self.peptides

    def getPeptideListOutputString(self):
        posList = []
        for nextPeptide in self.peptides.values():

            posList.append(nextPeptide.getOutputString())
        return '\n'.join(posList)
         
    def writeSequenceToFasta(self, outputDir):
        
        fastaFileName = os.path.join(outputDir, "%s.fasta" % self.modbaseSequenceId)
        fh = open(fastaFileName, 'w')
        fh.write(">%s\n" % self.modbaseSequenceId)

        fh.write("%s" % self.proteinSequence)
        fh.close()
        return fastaFileName

    def processDisopred(self, disopredReader, disopredRunner):
        """Run disopred if needed and read result file. Propogate results to peptides as features"""
        if self.hasErrors():
            return
        try:
            self.runDisopred(disopredReader, disopredRunner)
        except pcssErrors.DisopredException as pe:
            print "caught exception %s" % pe.msg
            for nextPeptide in self.peptides.values():
                print "next peptide %s process excpetion" % nextPeptide.startPosition
                nextPeptide.processException(pe)
            log.error(pe.msg)
            return
        self.populatePeptideDisopredValues()

    def runDisopred(self, disopredReader, disopredRunner):
        """Get Disopred calls for this protein and validate result file matches my sequence"""
        self.disopredProteinCalls = self.loadProteinSequenceFeatures(disopredReader, disopredRunner)
        self.disopredProteinCalls.checkSequenceMatch(self.proteinSequence)

    def populatePeptideDisopredValues(self):
        """Initialize disopred values for all my peptides"""
        for nextPeptide in self.peptides.values():
            try:
                nextPeptide.setDisopredResult(self.disopredProteinCalls)
            except pcssErrors.DisopredException as pe:
                nextPeptide.processException(pe)

    def processPsipred(self, psipredReader, psipredRunner):
        """Run psipred if needed and read result file. Propogate results to peptides as features"""
        if self.hasErrors():
            return
        try:
            self.runPsipred(psipredReader, psipredRunner)
        except pcssErrors.PsipredException as pe:
            for nextPeptide in self.peptides.values():
                nextPeptide.processException(pe)
            log.error(pe.msg)
            return
        self.populatePeptidePsipredValues()

    def runPsipred(self, psipredReader, psipredRunner):
        """Get Psipred calls for this protein and validate result file matches my sequence"""
        self.psipredProteinCalls = self.loadProteinSequenceFeatures(psipredReader, psipredRunner)
        self.psipredProteinCalls.checkSequenceMatch(self.proteinSequence)

    def populatePeptidePsipredValues(self):
        """Initialize psipred values for all my peptides"""
        for nextPeptide in self.peptides.values():
            try:
                nextPeptide.setPsipredResult(self.psipredProteinCalls)
            except pcssErrors.PsipredException as pe:
                nextPeptide.processException(pe)

    def loadProteinSequenceFeatures(self, sequenceReader, sequenceRunner):
        """Generic method for running a protein sequence feature algorithm and reading the result"""
        if (not sequenceReader.sfh.outputFileExists(self.modbaseSequenceId)):
            log.debug("Protein %s: %s result file does not exist; creating new..." % (self.modbaseSequenceId, 
                                                                                       sequenceReader.sfh.getName()))
            sequenceRunner.runSequenceFeature(self)
        else:
            log.debug("Protein %s: %s result file exists; reading..." % (self.modbaseSequenceId, 
                                                                         sequenceReader.sfh.getName()))
        sequenceFeatureResult = sequenceReader.readResult(self.modbaseSequenceId)

        return sequenceFeatureResult

    def addModels(self, modelTable):
        """Get models for my protein from model table object and add it to my internal dictionary"""
        if self.hasErrors():
            return
        
        modelSequence = modelTable.getPcssModelSequence(self.modbaseSequenceId)
        print "got model squence; %s" % modelSequence
        modelList = modelSequence.getModels()
        print "model list length: %s" % len(modelList)
        self._pcssModels = {}
        for model in modelList:
            model.calculateCoverage(self.getSequenceLength())
            self._pcssModels[model.getAttributeValue("model_id")] = model

    def processDssp(self):
        """Run DSSP on each of my peptide's best model"""
        if self.hasErrors():
            return
        for nextPeptide in self.peptides.values():
            if (self.hasModels()):
                nextPeptide.setBestModel(self.getRankedModels())
                try:
                    #nextPeptide.setTemplate()
                    nextPeptide.processDssp()
                except pcssErrors.StructureException as e:
                    nextPeptide.processException(e)
                    continue
            else:
                nextPeptide.addEmptyDsspFeatures()

    def validatePeptideSequences(self):
        for peptide in self.peptides.values():

            if (peptide.sequence != self.getSubsequence(peptide.startPosition, peptide.endPosition + 1)):
                raise pcssErrors.PcssGlobalException("Protein %s subsequence %s doesn't match peptide sequence %s starting at position %s" %
                                                     (self.modbaseSequenceId, self.getSubsequence(peptide.startPosition, peptide.endPosition + 1), 
                                                      peptide.sequence, peptide.startPosition))

    def getSubsequence(self, start, end):

        return str(self.proteinSequence[start:end])

    def getRankedModels(self):
        """Return a list of models ranked by the user-input best model criteria"""
        return sorted (self._pcssModels.values(), 
                       key=lambda model: model.getAttributeValue(self.pcssRunner.pcssConfig["best_model_attribute"]), reverse=True)
            

    def isEqual(self, otherProtein):
        for attName in otherProtein.proteinAttributes.keys():
            if (attName not in self.proteinAttributes):
                print "protein ie1"
                return False
            if (str(self.getAttributeOutputString(attName)) != str(otherProtein.getAttributeOutputString(attName))):
                print "protein ie2 attriute %s" % attName
                return False
        for otherPeptide in otherProtein.peptides.values():
            if (otherPeptide.startPosition not in self.peptides):
                print "protein ie3"
                return False
            myPeptide = self.peptides[otherPeptide.startPosition]
            if (not otherPeptide.isEqual(myPeptide)):
                print "protein ie4"
                return False
        return True
            

class PcssPeptide:
    
    """Class for one peptide; provides feature tracking and conversion methods"""

    def __init__(self, sequence, startPosition, endPosition, pcssRunner):
        self.sequence = sequence
        self.startPosition = startPosition
        self.endPosition = endPosition
        self.bestModel = None
        self.attributes = {}
        self.pcssRunner = pcssRunner
        self.addFeature(pcssFeatures.PeptideSequenceFeature(sequence))
        self.addStringAttribute("peptide_start", startPosition)
        self.addStringAttribute("peptide_end", endPosition)
        self.addFeature(pcssFeatures.PeptideErrorFeature())

    def getPeptideLength(self):
        return len(self.sequence)

    def getAttributeOutputString(self, attributeName):
        """Return the requested attribute as a string

        In contrast to PcssProteins, PcssPeptides have all attributes represented as PeptideFeature objects
        as they can be slightly more complex. This method returns the string value of that object as opposed
        to the feature object itself (which getAttribute() is used for)"""
        if (not self.hasAttribute(attributeName)):
            return None        
        return self.getAttribute(attributeName).getValueString()

    def getAttribute(self, attributeName):
        """Return the requested attribute object"""
        if (not self.hasAttribute(attributeName)):
            return None
        else:
            return self.attributes[attributeName]

    def hasAttribute(self, attributeName):
        return attributeName in self.attributes

    def getAttributeNameString(self):
        return ", ".join(self.attributes.keys())

    def addError(self, errorCode):
        self.getAttribute("peptide_errors").addError(errorCode)

    def addStringAttribute(self, name, value):
        """Quick way to add an attribute that is only a string value"""
        self.pcssRunner.pfa.validateAttribute(name)
        stringAttribute = pcssFeatures.StringAttribute(name, value)
        self.addFeature(stringAttribute)

    def getRange(self):
        """Return range as a two-element list"""
        return [self.startPosition, self.endPosition]

    def processException(self, e):
        """Store information from this exception internally
        
        Exceptions for peptides generally apply to a set of PeptideFeatures. To mark those as unavailable,
        the exception provides their attribute names. Here, set those values to be the error code from the exception"""
        affectedAttNames = e.getAffectedAttributeNames()
        for attName in affectedAttNames:
            self.addStringAttribute(attName, e.code)
        self.addError(e.code)

    def createDsspFeatures(self):
        """Get solvent accessibility and secondary structure as assessed by DSSP and create respective features for them"""
        dsspAccList = []
        dsspStructureList = []
        for i in range(len(self.sequence)):
            nextResidue = self.sequence[i]
            nextPosition = self.startPosition + i
            if (nextResidue != self.bestModel.getDsspResidueCode(nextPosition)):
                raise pcssErrors.DsspMismatchException("Dssp residue %s does not match peptide sequence %s at position %s" %
                                                       (self.bestModel.getDsspResidueCode(nextPosition), nextResidue, nextPosition))
            dsspAccList.append(self.bestModel.getRelativeSolventAcc(nextPosition))
            dsspStructureList.append(self.bestModel.getSecondaryStructure(nextPosition))
        dsspStructureFeature = pcssFeatures.DsspStructureFeature(dsspStructureList)
        dsspAccFeature = pcssFeatures.DsspAccFeature(dsspAccList)

        self.addFeature(dsspStructureFeature)
        self.addFeature(dsspAccFeature)

    def addEmptyDsspFeatures(self):
        self.addFeature(pcssFeatures.DsspStructureFeature())
        self.addFeature(pcssFeatures.DsspAccFeature())


    def processDssp(self):
        """Run DSSP for my best model and create features based on the results"""
        if (not self.bestModel):
            self.addEmptyDsspFeatures()
            return
        try:
            self.bestModel.loadDsspResults()
            self.createDsspFeatures()
        except pcssErrors.StructureException as e:
            e.setPeptide(self)
            e.setModel(self.bestModel)
            raise e

    def getOutputString(self):
        output =  "Start: %s end %s sequence %s\n" % (self.startPosition, self.endPosition, self.sequence)
        for feature in self.attributes.values():
            output += "%s\n" % feature.getOutputString()
        return output

    def setDisopredResult(self, disorderProteinCalls):
        """Use DisopredSequenceFeatureCallSet for full protein sequence to create disorder features for this peptide"""
        disorderStringList = []
        disorderScoreList = []
        for i in range(self.startPosition, self.endPosition + 1):
            try:
                nextCall = disorderProteinCalls.getSequenceFeatureCall(i+1)
                disorderStringList.append(nextCall.call)
                disorderScoreList.append(nextCall.score)
            except pcssErrors.DisopredPeptideNotFoundException as e:
                e.setPeptide(self)
                raise e
                                
        disorderStringFeature = pcssFeatures.DisorderStringFeature(disorderStringList)
        disorderScoreFeature = pcssFeatures.DisorderScoreFeature(disorderScoreList)

        self.addFeature(disorderStringFeature)
        self.addFeature(disorderScoreFeature)

    def setPsipredResult(self, psipredProteinCalls):
        """Use PsipredSequenceFeatureCallSet for full protein sequence to create disorder features for this peptide"""
        psipredStringList = []
        psipredScoreList = []
        for i in range(self.startPosition, self.endPosition + 1):
            try:
                nextCall = psipredProteinCalls.getSequenceFeatureCall(i+1)
                psipredStringList.append(nextCall.call)
                psipredScoreList.append(nextCall.score)
            except pcssErrors.PsipredPeptideNotFoundException as e:
                e.setPeptide(self)
                raise e
                                
        psipredStringFeature = pcssFeatures.PsipredStringFeature(psipredStringList)
        psipredScoreFeature = pcssFeatures.PsipredScoreFeature(psipredScoreList)

        self.addFeature(psipredStringFeature)
        self.addFeature(psipredScoreFeature)

    def addFeature(self, feature):
        self.pcssRunner.pfa.validateAttribute(feature.name)
        self.attributes[feature.name] = feature

    def setTemplate(self):
        if (not self.bestModel):
            return

    def setBestModel(self, modelList):
        for model in modelList:
            if (model.containsPeptide(self)):
                self.bestModel = model
                break

    def hasBestModel(self):
        return self.bestModel is not None

    def isEqual(self, otherPeptide):
        
        for attName in otherPeptide.attributes.keys():
            if (attName not in self.attributes and otherPeptide.getAttribute(attName).isInitialized()):
                print "peptide ie1"
                return False
            if (str(self.getAttributeOutputString(attName)) != str(otherPeptide.getAttributeOutputString(attName))):
                print "peptide start %s ie2 attribute %s" % (self.startPosition, attName)
                return False
        if (self.hasBestModel() != otherPeptide.hasBestModel()):
            print "peptide ie3"
            return False
        if (self.hasBestModel()):
            if (not otherPeptide.bestModel.isEqual(self.bestModel)):
                print "peptide ie4"
                return False
        return True

    def makeSvmFileLine(self):
        svmFileStringList = []
        featureNumber = 0
        svmHandler = pcssSvm.SvmFeatureHandler()
        self.svmHandler = svmHandler
        for featureName in self.pcssRunner.getSvmFeatureOrder():
            if (not self.hasAttribute(featureName)):
                raise pcssErrors.PcssGlobalException("Error: peptide tried to make svm feature for %s but does not have this feature" % featureName)
            
            if (not self.getAttribute(featureName).isInitialized() or pcssTools.isPeptideErrorValue(self.getAttributeOutputString(featureName))):
                svmHandler.processEmptyFeature(self.pcssRunner.getPeptideLength(), self.getAttribute(featureName))
                
            else:
                svmFileStringList.append(self.getAttribute(featureName).makeSvmFeature(svmHandler))
            svmHandler.finalizeFeature(self, self.getAttribute(featureName), self.pcssRunner.getPeptideLength())
        return " ".join(svmFileStringList)
                                     
