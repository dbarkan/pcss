import sys
import os
import re
import subprocess
from Bio import SeqIO
import configobj
import pcssTools
import pcssPeptide
import pcssErrors
import pcssIO
import pcssModels
import pcssFeatures
import pcssFeatureHandlers
import logging

#quick and dirty script to read annotation file and randomly assign positive and negative status
#to create defined input fasta file for testing


logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)


def makeModbaseIdsToSeqs(fastaFileName):
    fh = open(fastaFileName, 'r')
    fastaIterator = SeqIO.FastaIO.FastaIterator(fh)

    modbaseIdsToSeqs = {}
    for seqRecord in fastaIterator:
        [modbaseSeqId, uniprotId] = seqRecord.id.split('|')
        modbaseIdsToSeqs[modbaseSeqId] = str(seqRecord.seq)
    return modbaseIdsToSeqs

if (len(sys.argv) < 6):
    print "Usage: python makeDefinedInputFromAttributes.py <configFile> <configSpecFile> <inputFileName> <fastaFile> <outputFileName>"
    sys.exit()

configFile = sys.argv[1] #"testConfig/testPcssConfig.txt"
configSpecFile = sys.argv[2] #"testConfig/testConfigSpec.txt"
inputFileName = sys.argv[3]
fastaFile = sys.argv[4]
outputFileName = sys.argv[5]


try:
    pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)
    runner = pcssTools.PcssRunner(pcssConfig)

    reader = pcssIO.AnnotationFileReader(runner)
    reader.readAnnotationFile(inputFileName)

    modbaseIdsToSeqs = makeModbaseIdsToSeqs(fastaFile)

    proteins = reader.getProteins()
    peptideCount = 0
    outputFh = open(outputFileName, 'w')

    for protein in proteins:
        peptides = protein.peptides.values()
        peptideCodeList = []
        for peptide in peptides:
            status = ""
            if (peptideCount % 3 == 0):
                status = "positive"
            else:
                status = "negative"
            peptideCount += 1
            peptideCode = "_".join([str(peptide.startPosition), peptide.sequence, status])
            peptideCodeList.append(peptideCode)
        outputList = [protein.modbaseSequenceId, protein.uniprotId, "|".join(peptideCodeList)]
        outputFh.write(">%s\n"  % "|".join(outputList))
        outputFh.write("%s\n" % modbaseIdsToSeqs[protein.modbaseSequenceId])
except pcssErrors.PcssGlobalException as e:
    print e.msg
    
