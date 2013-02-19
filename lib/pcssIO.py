import sys
import os
from Bio import SeqIO
import itertools
import pcssModels
import pcssTools
import pcssCluster
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

        
        self.pcssRunner = pcssRunner

    def readInputFile(self, proteinFastaFile):
        """Read input fasta file, parse headers to create PcssProteins and parse sequence to create PcssPeptides"""
        print "DEFINED: Reading input %s" % proteinFastaFile
        
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
        print "making all peptides from seq id %s" % modbaseSeqId
        for col in cols[2:len(cols)]:
            nextPeptide = self.makePeptideFromCode(col)
            if (nextPeptide is not None):
                peptideList.append(nextPeptide)
        
        if (len(peptideList) < 1):
            raise pcssErrors.PcssGlobalException("Protein %s has no peptides" % modbaseSeqId)
        seq.setPeptides(peptideList)
        return seq

    def makePeptideFromCode(self, peptideCode):
        print "making peptide from %s" % peptideCode
        [peptideStart, peptideSequence, status] = peptideCode.split('_')
        if (peptideSequence == self.pcssRunner.internalConfig["keyword_peptide_sequence_mismatch"]):
            return None
        status = self.pcssRunner.validatePeptideCodeStatus(status, peptideCode)
        peptideStart = int(peptideStart)
        peptide = pcssPeptide.PcssPeptide(peptideSequence, peptideStart, peptideStart + len(peptideSequence) -1, self.pcssRunner)
        peptide.addStringAttribute("status", status)
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
                self.validateColumnLine(annotationFile, line)

                continue
            pcssProtein = self.getProteinFromLine(line)
            if (not(pcssProtein.hasErrors())):
                cols = line.split('\t')
                for attribute in sortedAttributes:
                    attribute.setValueFromFile(self.getValueForAttributeName(attribute.name, cols), 
                                               pcssProtein, 
                                               int(self.getValueForAttributeName("peptide_start", cols)))
        
        if (len(self.proteins) == 0):
            raise pcssErrors.PcssGlobalException("Did not read any proteins from annotation file")

    def readProteinSequences(self, fastaFileName):
        fh = open(fastaFileName, 'r')
        fastaIterator = SeqIO.FastaIO.FastaIterator(fh)
        for seqRecord in fastaIterator:
            [modbaseId, uniprotId] = seqRecord.id.split('|')
            if (modbaseId in self.proteins):
                protein = self.proteins[modbaseId]
                protein.setProteinSequence(seqRecord.seq)

        for protein in self.proteins.values():
            if (protein.proteinSequence is None):
                raise pcssErrors.PcssGlobalException("Protein %s has no sequence set" % protein.modbaseSequenceId)
    def validateColumnLine(self, annotationFile, line):
        sortedAttributes = self.pcssRunner.pfa.getColumnSortedInputAttributes()
        firstAttribute = sortedAttributes[0]
        if (not line.startswith(firstAttribute.niceName)):
            raise pcssErrors.PcssGlobalException("Error: read annotation file %s\n. Expected first row to be column header "
                                                 "(starting with %s) but didn't find it; instead got\n%s" % (annotationFile, 
                                                                                                             firstAttribute.niceName, line))
        columnNames = line.split('\t')
        sortedAttributeNames = []
        for i in sortedAttributes:
            sortedAttributeNames.append(i.niceName)
            if (i.niceName not in columnNames):
                raise pcssErrors.PcssGlobalException("Error: read annotation file %s\n. Expected input attribute %s but did not find it"
                                                     % (annotationFile, i.niceName))
        for i in columnNames:
            if i not in sortedAttributeNames:
                raise pcssErrors.PcssGlobalException("Error: read annotation file %s\n. Read column header %s that wasn't specified in attributes file"
                                                     % (annotationFile, i))
        
    def getProteinFromLine(self, line):
        cols = line.split('\t')
        

        protein = self.getProtein(cols)
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
         return cols[self.pcssRunner.pfa.getAttribute(attributeName).inputOrder]
            
    def getProtein(self, cols):
        sequenceId = self.getValueForAttributeName("seq_id", cols)

        if (sequenceId not in self.proteins):
            protein = pcssPeptide.PcssProtein(sequenceId, self.pcssRunner)
            self.proteins[sequenceId] = protein
            uniprotId = self.getValueForAttributeName("uniprot_id", cols)
            protein.setStringAttribute("uniprot_id", uniprotId)
            protein.uniprotId = self.getValueForAttributeName("uniprot_id", cols)

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
    def __init__(self, name, attributeType, optional, niceName, featureClass, io):
        self.name = name
        self.attributeType = attributeType
        if (optional == "True"):
            self.outputOptional = True
        else:
            self.outputOptional = False
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

    def setInputOrder(self, inputOrder):
        self.inputOrder = inputOrder

    def setOutputOrder(self, outputOrder):
        self.outputOrder = outputOrder

    def isInputAttribute(self):
        return self.input

    def isOutputAttribute(self):
        return self.output

    def getProteinValue(self, protein):
        """Get value of my attribute from input protein as a string"""

        attributeValue = protein.getAttributeOutputString(self.name)

        if (attributeValue is None): 
            if (self.outputOptional is False):
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
                if (self.outputOptional is False):
                    raise pcssErrors.PcssGlobalException("Peptide %s never set mandatory attribute %s" % (peptide.startPosition, self.name))
                else:
                    return ""
            return attributeValue

    def setValueFromFile(self, fileValue, protein, peptideStartPosition):

        if (self.attributeType == "protein"):
            protein.setStringAttribute(self.name, fileValue) #currently all attributes are string attributes
        elif(self.attributeType == "peptide"):
            peptide = protein.peptides[peptideStartPosition]
            if (self.isError(fileValue)):
                peptide.addStringAttribute(self.name, fileValue)
            else:
                classObject = self.getFeatureClassObject()
                if (self.featureClass == "StringAttribute"):
                    classObject.initFromFileValue(self.name, fileValue)
                else:
                    classObject.initFromFileValue(fileValue)
                peptide.addFeature(classObject)
        else:
            if (fileValue != ""):
                protein.peptides[peptideStartPosition].bestModel.setAttribute(self.name, fileValue)

    def isError(self, attributeValue):
        return attributeValue.startswith("peptide_")
    

    def getFeatureClassObject(self):
        fullClassName = "pcssFeatures.%s" % self.featureClass
        className = eval(fullClassName)
        classObject = className()
        return classObject

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
        inputCounter = 0
        outputCounter = 0
        for line in lines:

            [name, attributeType, optional, niceName, featureClass, io] = line.split('\t')
            att = PcssFileAttribute(name, attributeType, optional, niceName, featureClass, io)
            if (att.isInputAttribute()):
                att.setInputOrder(inputCounter)
                inputCounter += 1
            if (att.isOutputAttribute()):
                att.setOutputOrder(outputCounter)
                outputCounter += 1
            self.setFileAttribute(att)

    def setFileAttribute(self, attribute):
        self._attributes[attribute.name] = attribute

    def getAttribute(self, name):
        return self._attributes[name]

    def setAllOptional(self):
        """Set all attributes to be optional regardless of what was in the file; useful for testing"""
        for att in self._attributes.values():
            att.outputOptional = True

    def getAttributeListByType(self, targetType):
        attributeList = []
        for attribute in self._attributes.values():
            if (attribute.attributeType == targetType):
                attributeList.append(attribute)
        return attributeList


    def getColumnSortedInputAttributes(self):
        """Get all input attributes objects sorted by their output order"""
        inputAttributes = []
        for attribute in self._attributes.values():
            if (attribute.isInputAttribute()):
                inputAttributes.append(attribute)
        return sorted (inputAttributes, key=lambda attribute: attribute.inputOrder)

    def getColumnSortedOutputAttributes(self):
        """Get all input attributes objects sorted by their output order"""
        outputAttributes = []
        for attribute in self._attributes.values():
            if (attribute.isOutputAttribute()):
                outputAttributes.append(attribute)
        return sorted (outputAttributes, key=lambda attribute: attribute.outputOrder)
        

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
        self.outputFh = open(pcssRunner.pdh.getFullOutputFile(pcssRunner.internalConfig["annotation_output_file"]), 'w')

    def writeGlobalException(self):
        print "global exception"

    def writeAllOutput(self, proteins):
        """Write output file; write column headers and write one line for each peptide in the protein set."""
        self.outputFh.write("%s\n" % self.pcssRunner.pfa.getOutputColumnHeaderString())
        for protein in proteins:
            self.writeProteinOutputLines(protein)
        self.outputFh.close()

    def writeProteinOutputLines(self, protein):
        if (protein.hasErrors()):
            self.outputFh.write(self.makeProteinErrorLine(protein))
        else:
            peptides = protein.peptides.values()
            for peptide in peptides:
                self.outputFh.write(self.makePeptideOutputLine(protein, peptide))

    def makeProteinErrorLine(self, protein):
        outputList = self.makeProteinErrorOutputList(protein)        
        outputString = "%s\n" % '\t'.join(str(x) for x in outputList)
        return outputString

    def makePeptideOutputLine(self, protein, peptide):
        """Return a string representing the output line for this peptide"""
        outputList = self.makePeptideOutputList(protein, peptide)
        outputString = "%s\n" % '\t'.join(str(x) for x in outputList)
        return outputString

    #def makeOutputList(self, protein, peptide):
    #    """Return a list of protein and peptide attributes, ordered by specified file attribute order.

    #    If protein has errors affecting the whole protein (e.g., doesn't have any peptides) then return a list
    #    containing only the protein sequence and error code, and empty elements where other attributes would 
    #    normally be"""
    #    return self.makeProteinOutputList(protein, peptide)

    def makePeptideOutputList(self, protein, peptide):
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
            if (attribute.name == "seq_id" or attribute.name == "protein_errors" or attribute.name == "uniprot_id"):
                #sort of a hack to include uniprot_id; don't need it in output list but need it to get the tests to pass
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



