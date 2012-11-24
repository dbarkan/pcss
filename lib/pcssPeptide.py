import sys
import os
from Bio import SeqIO

class PcssProtein:

    """Class for one protein sequence, containing a number of peptides on which features are assessed"""
    
    def __init__(self, modbaseSequenceId):
        self.modbaseSequenceId = modbaseSequenceId
        self.peptides = {}

class PcssPeptide:
    
    """Class for one peptide; provides feature tracking and conversion methods"""

    def __init__(self, sequence, startPosition, endPosition):
        self.sequence = sequence
        self.startPosition = startPosition
        self.endPosition = endPosition

        self.features = {}

class ScanPeptideImporter:
    def readPeptides(self, peptideFile):
        fh = open(peptideFile, 'r')
        fastaIterator = SeqIO.FastaIO.FastaIterator(fh)
        for seqRecord in fastaIterator:
            print seqRecord.id
            print seqRecord.seq

#class DefinedPeptideImporter:
    
    

    
