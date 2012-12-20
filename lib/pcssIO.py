import sys
import os
from Bio import SeqIO
import itertools
import pcssModels
import pcssTools

import logging
import pcssPeptide
import pcssErrors
import pcssFeatures
import pcssFeatureHandlers
log = logging.getLogger("pcssPeptide")

class ScanPeptideImporter:

    """Class to read a fasta sequence and parse peptides from the actual sequence according to defined rules"""

    def __init__(self, pcssRunner):

        self.rules = ParsingRules(pcssRunner.pcssConfig['rules_file'])
        self.peptideLength = pcssRunner.pcssConfig.as_int('peptide_length')
        self.pcssRunner = pcssRunner

    def readInputFile(self, proteinFastaFile):
        """Read input fasta file, parse headers to create PcssProteins and parse sequence to create PcssPeptides"""
        fh = open(proteinFastaFile, 'r')
        fastaIterator = SeqIO.FastaIO.FastaIterator(fh)
        pcssProteins = []
        
        for seqRecord in fastaIterator:
            
            pcssProtein = self.parseFastaHeader(seqRecord)
            pcssProtein.setProteinSequence(seqRecord.seq)
            pcssPeptideList = self.parseFastaSequence(seqRecord)
            pcssProtein.setPeptides(pcssPeptideList)

            pcssProteins.append(pcssProtein)
        
        log.info("ScanPeptideImporter: read %s proteins from input file" % len(pcssProteins))
        return pcssProteins
        
    def parseFastaHeader(self, seqRecord):
        """Return PcssProtein by reading fasta header. Header is of format <modbaseId|uniprotId>"""
        [modbaseId, uniprotId] = seqRecord.id.split('|')
        seq = pcssPeptide.PcssProtein(modbaseId, self.pcssRunner)
        seq.setUniprotId(uniprotId)
        return seq

    def parseFastaSequence(self, seqRecord):
        """Use user-defined rules to read protein sequence, parse peptides, and return a list of those that conform to rules"""
        seqLength = len(seqRecord.seq)
        pcssPeptideList = []
        for i in range(0, seqLength - self.peptideLength + 1):
            nextPeptide = str(seqRecord.seq[i:i+self.peptideLength])
            if (self.rules.isValidPeptide(nextPeptide)):
                pcssPeptideList.append(pcssPeptide.PcssPeptide(nextPeptide, i, i + self.peptideLength - 1, self.pcssRunner))
        return pcssPeptideList


class DefinedPeptideImporter:

    """Class to read a fasta file where peptides to process are defined in the header for each sequence"""

    def __init__(self, pcssRunner):

        self.rules = ParsingRules(pcssRunner.pcssConfig['rules_file'])
        self.peptideLength = pcssRunner.pcssConfig.as_int('peptide_length')
        self.pcssRunner = pcssRunner

    def readInputFile(self, proteinFastaFile):
        """Read input fasta file, parse headers to create PcssProteins and parse sequence to create PcssPeptides"""
        fh = open(proteinFastaFile, 'r')
        fastaIterator = SeqIO.FastaIO.FastaIterator(fh)
        pcssProteins = []

        for seqRecord in fastaIterator:
            
            pcssProtein = self.parseFastaHeader(seqRecord)
            pcssProtein.setProteinSequence(seqRecord.seq)
            pcssProtein.validatePeptideSequences()
            pcssProteins.append(pcssProtein)
            
        
        log.info("Read %s proteins from input file" % len(pcssProteins))
        return pcssProteins
        
    def parseFastaHeader(self, seqRecord):
        cols = seqRecord.id.split('|')
        modbaseSeqId = cols[0]
        uniprotId = cols[1]
        seq = pcssPeptide.PcssProtein(modbaseSeqId, self.pcssRunner)
        seq.setUniprotId(uniprotId)
    
        peptideList = []
        for col in cols[2:-1]:
            peptideList.append(self.makePeptideFromCode(col))
        
        seq.setPeptides(peptideList)
        return seq

    def makePeptideFromCode(self, peptideCode):
        [peptideStart, peptideSequence, classification] = peptideCode.split('_')
        peptideStart = int(peptideStart)
        peptide = pcssPeptide.PcssPeptide(peptideSequence, peptideStart, peptideStart + len(peptideSequence) -1, self.pcssRunner)
        return peptide

class AnnotationFileReader:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.proteins = {}

    def readAnnotationFile(self, annotationFile):
        if (not os.path.exists(annotationFile)):
            raise pcssErrors.PcssGlobalException("Error: annotation file reader did not find expected annotation file\n%s" % annotationFile)
        reader = pcssTools.PcssFileReader(annotationFile)
        lines = reader.getLines()
        sortedAttributes = self.pcssRunner.pfa.getColumnSortedInputAttributes()

        for (i, line) in enumerate(lines):
            if (i == 0):
                firstAttribute = sortedAttributes[0]
                if (not line.startswith(firstAttribute.niceName)):
                    raise pcssErrors.PcssGlobalException("Error: read annotation file %s\n. Expected first row to be column header "
                                                         "(starting with %s) but didn't find it; instead got\n%s" % (annotationFile, 
                                                                                                                     firstAttribute.niceName, line))
                continue
            pcssProtein = self.getProteinFromLine(line)
            if (not(pcssProtein.hasErrors())):
                cols = line.split('\t')
                for attribute in sortedAttributes:
                    attribute.setValueFromFile(cols[attribute.order], pcssProtein, int(self.getValueForAttributeName("peptide_start", cols)))
            
    def getProteinFromLine(self, line):
        cols = line.split('\t')
        
        sequenceId = self.getValueForAttributeName("seq_id", cols)
        protein = self.getProtein(sequenceId)
        protein.setStringAttribute("protein_errors", self.getValueForAttributeName("protein_errors", cols))
        
        if (not protein.hasErrors()):
            peptide = self.getPeptideFromLine(cols)
            model = self.getModelFromLine(cols)
            protein.setPeptide(peptide)
            peptide.bestModel = model
                        
        return protein

    def getModelFromLine(self, cols):
        modelId = self.getValueForAttributeName("model_id", cols)
        model = None
        if (modelId != ""):
            model = pcssModels.PcssModel(self.pcssRunner)
        return model

    def getPeptideFromLine(self, cols):
        
        return pcssPeptide.PcssPeptide(self.getValueForAttributeName("peptide_sequence", cols),
                                       int(self.getValueForAttributeName("peptide_start", cols)),
                                       int(self.getValueForAttributeName("peptide_end", cols)),
                                       self.pcssRunner)
    
                             
    def getValueForAttributeName(self, attributeName, cols):
         return cols[self.pcssRunner.pfa.getAttribute(attributeName).order]
            
    def getProtein(self, sequenceId):
        if (sequenceId not in self.proteins):
            self.proteins[sequenceId] = pcssPeptide.PcssProtein(sequenceId, self.pcssRunner)         
        return self.proteins[sequenceId]


    def getProteins(self):
        return self.proteins.values()


class ParsingRules:
    def __init__(self, rulesFile):
        """Read rules file and store user-defined rules for parsing peptides

        @paramRulesFile: full path of text file specifying user-defined rules for what peptides are
        of interest. Each line represents one position in the peptide. For each line, the first entry 
        is the position index and subsequent entries are one letter codes of amino acids not allowed 
        at that position. For example, the line "3 A C D" means that at position 3, a peptide will not
        have A, C, or D present. 
        """
        self._rules = {}

        reader = pcssTools.PcssFileReader(rulesFile)
        lines = reader.getLines()
        for line in lines:
            self.addRule(line)

    def addRule(self, ruleLine):
        """Add to self._rules by reading next file line

        @param ruleLine: space separated columns from text, formatted as in class description
        """
        cols = ruleLine.split(' ')
        positionNumber = int(cols[0])
        self._rules[positionNumber] = {}
        for i in range(1, len(cols)):
            self._rules[positionNumber][cols[i].upper()] = 1
                                        
    def isValidPeptide(self, sequence):
        """Return True if passed sequence conforms to rules"""
        for position, disallowedAAs in self._rules.iteritems():
            nextAA = sequence[position - 1].upper()
            if nextAA in disallowedAAs:
                return False
        return True

class PcssFileAttribute:

    """Class specifying the properties of proteins and peptides that will be written to result file

    This class manages the order that attributes are written to the file, whether they are mandatory
    or optional, and other properties."""
    def __init__(self, order, name, attributeType, optional, niceName, featureClass, io):
        self.name = name
        self.order = order
        self.attributeType = attributeType
        if (optional == "True"):
            self.optional = True
        else:
            self.optional = False
        self.niceName = niceName
        self.featureClass = featureClass
        ioCols = io.split(',')
        self.input = False
        self.output = False
        for col in ioCols:
            if (col == "input"):
                self.input = True
            if (col == "output"):
                self.output = True

    def isInputAttribute(self):
        return self.input

    def isOutputAttribute(self):
        return self.output

    def getProteinValue(self, protein):
        """Get value of my attribute from input protein as a string"""

        attributeValue = protein.getAttributeOutputString(self.name)

        if (attributeValue is None): 
            if (self.optional is False):
                raise pcssErrors.PcssGlobalException("Protein %s never set mandatory attribute %s" % (protein.modbaseSequenceId, self.name))
            else:
                return ""
        return attributeValue

    def getPeptideValue(self, peptide):
        """Get value of my attribute from input peptide as a string"""
        if (self.attributeType == 'model'):

            if (peptide.bestModel is None):
                return ""
            else:
                return peptide.bestModel.getAttributeValue(self.name)
        else:
            attributeValue = peptide.getAttributeOutputString(self.name)
            if (attributeValue is None): 
                if (self.optional is False):
                    raise pcssErrors.PcssGlobalException("Peptide %s never set mandatory attribute %s" % (peptide.startPosition, self.name))
                else:
                    return ""
            return attributeValue

    def setValueFromFile(self, fileValue, protein, peptideStartPosition):
        #if (fileValue != ""):
        if (self.attributeType == "protein"):
            protein.setStringAttribute(self.name, fileValue) #currently all attributes are string attributes
        elif(self.attributeType == "peptide"):
            peptide = protein.peptides[peptideStartPosition]
            fullClassName = "pcssFeatures.%s" % self.featureClass
            className = eval(fullClassName)

            classObject = className()
            if (self.featureClass == "StringAttribute"):
                classObject.initFromFileValue(self.name, fileValue)
            else:

                classObject.initFromFileValue(fileValue)
                peptide.addFeature(classObject)
        else:
            if (fileValue != ""):
                protein.peptides[peptideStartPosition].bestModel.setAttribute(self.name, fileValue)


class PcssFileAttributes:

    """Class for managing a set of file attributes"""

    def __init__(self, pcssConfig):
        self.pcssConfig = pcssConfig
        self.readAttributeFile()

    def readAttributeFile(self):
        """Read attribute table from input file. The order in which attributes are listed in the file are their output order"""
        reader = pcssTools.PcssFileReader(self.pcssConfig["attribute_file_name"])
        lines = reader.getLines()
        self._attributes = {}
        for (i, line) in enumerate(lines):

            [name, attributeType, optional, niceName, featureClass, io] = line.split('\t')
            att = PcssFileAttribute(i, name, attributeType, optional, niceName, featureClass, io)
            self.setFileAttribute(att)

    def setFileAttribute(self, attribute):
        self._attributes[attribute.name] = attribute

    def getAttribute(self, name):
        return self._attributes[name]

    def setAllOptional(self):
        """Set all attributes to be optional regardless of what was in the file; useful for testing"""
        for att in self._attributes.values():
            att.optional = True

    def getAttributeListByType(self, targetType):
        attributeList = []
        for attribute in self._attributes.values():
            if (attribute.attributeType == targetType):
                attributeList.append(attribute)
        return attributeList


    def getColumnSortedInputAttributes(self):
        """Get all input attributes objects sorted by their output order"""
        finalAttributes = []
        for att in self.getColumnSortedAttributes():
            if (att.isInputAttribute()):
                finalAttributes.append(att)
        return finalAttributes

    def getColumnSortedOutputAttributes(self):
        """Get all input attributes objects sorted by their output order"""
        finalAttributes = []
        for att in self.getColumnSortedAttributes():
            if (att.isOutputAttribute()):
                
                finalAttributes.append(att)
        return finalAttributes

    def getColumnSortedAttributes(self):
        return sorted (self._attributes.values(), key=lambda attribute: attribute.order)

    def validateAttribute(self, attributeName):
        """Check to see if this is an attribute originally specified in the input attribute table file"""
        if (not attributeName in self._attributes):
            raise pcssErrors.PcssGlobalException("Error: attempted to set attribute %s which is not a valid pfa attribute" % attributeName)
        
    def getOutputColumnHeaderString(self):
        return '\t'.join(x.niceName for x in self.getColumnSortedOutputAttributes())
        
class AnnotationFileWriter:
    
    """Class to write all protein and peptide output to a file"""

    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.outputFh = open(pcssRunner.pdh.getFullOutputFile(pcssRunner.pcssConfig["annotation_output_file"]), 'w')
            
    def writeGlobalException(self):
        print "global exception"

    def writeAllOutput(self, proteins):
        """Write output file; write column headers and write one line for each peptide in the protein set."""
        self.outputFh.write("%s\n" % self.pcssRunner.pfa.getOutputColumnHeaderString())
        for protein in proteins:
            peptides = protein.peptides.values()
            for peptide in peptides:
                self.writeOutputLine(protein, peptide)
        self.outputFh.close()

    def makeOutputLine(self, protein, peptide):
        """Return a string representing the output line for this peptide"""
        outputList = self.makeOutputList(protein, peptide)
        outputString = "%s\n" % '\t'.join(str(x) for x in outputList)
        return outputString

    def writeOutputLine(self, protein, peptide):
        """Write the string for this peptide to disk"""
        self.outputFh.write(self.makeOutputLine(protein, peptide))
        

    def makeOutputList(self, protein, peptide):
        """Return a list of protein and peptide attributes, ordered by specified file attribute order.

        If protein has errors affecting the whole protein (e.g., doesn't have any peptides) then return a list
        containing only the protein sequence and error code, and empty elements where other attributes would 
        normally be"""
        if (protein.hasErrors()):
            return self.makeProteinErrorOutputList(protein)
        else:
            return self.makeProteinOutputList(protein, peptide)

    def makeProteinOutputList(self, protein, peptide):
        """Return a list of protein and peptide attributes, ordered by specified file attribute order. Assumes no ProteinErrors"""
        attributeList = self.pcssRunner.pfa.getColumnSortedOutputAttributes()
        outputList = []
        for attribute in attributeList:
            if (attribute.attributeType == "protein"):
                outputList.append(attribute.getProteinValue(protein))
            else:
                outputList.append(attribute.getPeptideValue(peptide))
        return outputList

    def makeProteinErrorOutputList(self, protein):
        """Return a list with protein sequence and protein error code, with empty elements where other attributes would normally be."""
        attributeList = self.pcssRunner.pfa.getColumnSortedOutputAttributes()
        outputList = []
        for attribute in attributeList:
            if (attribute.name == "sequence_id" or attribute.name == "protein_errors"):
                outputList.append(attribute.getProteinValue(protein))
            else:
                outputList.append('')
        return outputList
                
class FinishPcss:
    def __init__(self):
        e = None # works?
        try:
            doStuff() #exceptions all handled within call
        except pcssErrors.PcssGlobalException as e:
            print e.msg
        finally:
            if e is not None:
                applicationFile.writeGlobalError
            else:
                applciationFile.output()
        
        

    #have one attribute for errors -- maybe global

    #we want to just get feature objects and for each write it out
    #can have an output line object in which feature objects are stored
    #also handle no peptides parsed

    #output line writes itself according to the order in the column file, can also write somethign different if there is an error
    #Errors:
    # Global error: output file writes the error code to first line and the message to the rest (no columns written)
    # full protein error: write protein sequence id and then the error code (should only be one and then no peptide columns so should be fine)
    # feature exception: write the error code to the peptide column, add to list in error column
    # peptides and proteins should save the exception as errors. maybe even features. calling code should die if not have peptide / protein set to make sure always have one



