import sys
import os
from Bio import SeqIO
import itertools
import pcssIO
import pcssTools
import logging
import pcssErrors
import pcssFeatures
import pcssFeatureHandlers
import collections

class BenchmarkResults:
    
    def __init__(self):
        self._results = []
        self.ScoreTuple = collections.namedtuple('svmScore', ['fpr', 'tpr', 'score'])

    def readBenchmarkFile(self, fileName):
        reader = pcssTools.PcssFileReader(fileName)
        lines = reader.getLines()
        for line in lines:
            cols = line.split()
            fpr = float(cols[0])
            tpr = float(cols[1])
            score = float(cols[2])
            self.validateScore(score)
            st = self.ScoreTuple(fpr, tpr, score)
            self._results.append(st)

    def validateScore(self, score):
        if (len(self._results) > 0):
            if (score > self._results[-1].score):
                raise pcssErrors.PcssGlobalException("Benchmark file error: got score %s that was larger than score from previous line %s" % 
                                                     (score, self._results[-1].score))
    def getClosestScoreTuple(self, score):
        for st in self._results:
            if (float(score) >  st.score):
                return st


class ClassifySvm:
    def __init__(self, pcssRunner, modelFile, benchmarkResultsFile):
        self.pcssRunner = pcssRunner
        self.modelFile = modelFile
        self.br = BenchmarkResults()
        self.br.readBenchmarkFile(benchmarkResultsFile)

    def setProteins(self, proteins):
        self.proteins = proteins
        self.peptides = []
        for protein in proteins:
            peptides = protein.peptides.values()
            for peptide in peptides:
                self.peptides.append(peptide)
        
    def writeClassificationFile(self):
        classificationFileName = self.getSvmInputFile()
        classificationFh = open(classificationFileName, 'w')
        for peptide in self.peptides:
            nextLine = peptide.makeSvmFileLine()
            classificationFh.write("%s %s\n" % ("0", nextLine))
        classificationFh.close()

    def classifySvm(self):
        svmCommandName = self.pcssRunner.internalConfig['svm_classify_command']
        classificationFileName = self.getSvmInputFile()

        if (not os.path.exists(classificationFileName)):
            raise pcssErrors.PcssGlobalException("Did not find test set file in expected location -- searched for\n%s" % classificationFileName)

        scoreFileName = self.getClassifyOutputFile()

        svmOutput = self.pcssRunner.pdh.runSubprocess([svmCommandName, classificationFileName, self.modelFile, scoreFileName])
        

    def readResultFile(self):
        resultFile = self.getClassifyOutputFile()

        if (not os.path.exists(resultFile)):
            raise pcssErrors.PcssGlobalException("Classify SVM could not read result file %s; \n"
                                                 "check to make sure svm_classify completed as suggested" % resultFile)
        reader = pcssTools.PcssFileReader(self.getClassifyOutputFile())
        lines = reader.getLines()
        if (len(lines) != len(self.peptides)):
            raise pcssErrors.PcssGlobalException("Result file has a different number of results (%s) than I have peptides (%s)" % 
                                                 (len(lines), len(self.peptides)))
        for (i, peptide) in enumerate(self.peptides):
            score = float(lines[i])
            st = self.br.getClosestScoreTuple(score)
        #peptide set feature for score, positive, negative

class ApplicationSvm(ClassifySvm):

    def getSvmInputFile(self):
        return self.pcssRunner.pdh.getSvmApplicationSetFile()
    
    def getClassifyOutputFile(self):
        return self.pcssRunner.pdh.getSvmApplicationOutputFile()

class SvmFeatureHandler:
    def __init__(self):
        self.featureNumber = 1

    def processEmptyFeature(self, peptide, feature):
        
        self.featureNumber += feature.getEmptyFeatureOffset(peptide)
        
    def getFeatureNumber(self):
        return self.featureNumber

    def processFeature(self, increment):
        self.featureNumber += increment
