import sys
import os
from Bio import SeqIO
import itertools
import math
import pcssIO
import pcssTools
import random
import logging
import pcssErrors
import pcssFeatures
import pcssFeatureHandlers
import collections

class BenchmarkResults:
    
    def __init__(self):
        self._results = []
        self.ScoreTuple = collections.namedtuple('svmScore', ['fpr', 'tpr', 'score'])

    def checkBoundaryLines(self, line, value):
        value = str(value)
        line = line.rstrip()
        cols = line.split('\t')
        if (cols[0] == value and cols[1] == value and len(cols) == 2):
            return True
        return False

    def readBenchmarkFile(self, fileName):
        reader = pcssTools.PcssFileReader(fileName)
        lines = reader.getLines()
        firstLine = lines[0]
        lastLine = lines[-1]
        if (not self.checkBoundaryLines(firstLine, 0)):
            raise pcssErrors.PcssGlobalException("Expected benchmark file %s to have first line of 0\t0")
        if (not self.checkBoundaryLines(lastLine, 1)):
            raise pcssErrors.PcssGlobalException("Expected benchmark file %s to have last line of 1\t1")
        
        for line in lines[1:len(lines) - 2]:
            
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

class SvmBenchmarker:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.trainingSvm = TrainingSvm(pcssRunner)
        self.testSvm = TestSvm(pcssRunner)  
        
    def getAllPeptides(self, proteins):
        return pcssTools.getAllPeptides(proteins, False)


    def createTrainingAndTestSets(self, peptides):
        
        print "got %s peptides" % len(peptides)
        self.benchmarkHandler = TrainingBenchmarkHandler(self.pcssRunner, peptides)
        self.benchmarkHandler.makeTrainingAndTestSets()


        print "%s positive training set, %s positive test set, %s negative training set, %s negative test set" % (len(self.benchmarkHandler.positiveTrainingSet),
                                                                                                                  len(self.benchmarkHandler.positiveTestSet),
                                                                                                                  len(self.benchmarkHandler.negativeTrainingSet),
                                                                                                                  len(self.benchmarkHandler.negativeTestSet))
        self.fullTrainingSet = self.benchmarkHandler.positiveTrainingSet + self.benchmarkHandler.negativeTrainingSet
        self.fullTestSet = self.benchmarkHandler.positiveTestSet + self.benchmarkHandler.negativeTestSet
        
        self.trainingSvm.setPeptides(self.fullTrainingSet)
        self.trainingSvm.writeTrainingSetFile()
        self.testSvm.setPeptides(self.fullTestSet)
        self.testSvm.writeClassificationFile()

    

    
    def trainAndApplyModel(self):
        
        self.trainingSvm.trainModel()
        
        self.testSvm.classifySvm()

    def readBenchmarkResults(self):
        self.testSvm.readResultFile()
        
        pstList = self.testSvm.getBenchmarkScoreTupleList()
        for i in range(len(self.fullTestSet)):
            nextTuple = pstList.getBenchmarkTuple(i)
            print "next peptide %s score %s fpr %s" % (nextTuple.peptide.startPosition, nextTuple.score, nextTuple.fpr)
        #self.pstTracker.addPstList(pstList)
        #peptide score tuple list:
        #should be able to tell me what the fpr, tpr is at each point
        #should be able to tell the critical point is
        #if I have a list, get the average of fpr, tpr and standard deviation
        #SinglePeptideScoreTupleList
        #get the fpr, tpr at a certain index

    
        
        
class TrainingSvm:
    def __init__(self, runner):
        self.runner = runner

    def setPeptides(self, peptides):
        for peptide in peptides:
            self.runner.validatePeptideTrainingStatus(peptide.getAttributeOutputString("status"))
        self.peptides = peptides

    def getStatusCode(self, peptide):
        
        if (peptide.getAttributeOutputString("status") == self.runner.getPositiveKeyword()):
            return 1
        else:
            return -1
    
    def writeTrainingSetFile(self):
        trainingSetFileName  = self.runner.pdh.getSvmTrainingSetFile()
        trainingSetFh = open(trainingSetFileName, 'w')
        for peptide in self.peptides:
            nextLine = peptide.makeSvmFileLine()
            statusCode = self.getStatusCode(peptide)
            trainingSetFh.write("%s %s\n" % (statusCode, nextLine))
        trainingSetFh.close()

    def trainModel(self):
        svmCommandName = self.runner.internalConfig['svm_train_command']
        trainingSetFileName =  self.runner.pdh.getSvmTrainingSetFile()
        if (not os.path.exists(trainingSetFileName)):
            raise pcssErrors.PcssGlobalException("Did not find training set input file in expected location -- searched for\n%s" % trainingSetFileName)
        modelFileName = self.runner.pdh.getSvmNewModelFile()
        gammaFlag = self.runner.pcssConfig["svm_training_gamma"]
        cFlag = self.runner.pcssConfig["svm_training_c"]
        
        #SPLIT FLAGS
        svmOutput = self.runner.pdh.runSubprocess([svmCommandName, "-g", gammaFlag, "-c", cFlag, trainingSetFileName, modelFileName])
        
    

class TrainingBenchmarkHandler:
    def __init__(self, pcssRunner, peptides):

        self.positivePeptides = self.makePeptideSet(peptides, pcssRunner.getPositiveKeyword())
        self.negativePeptides = self.makePeptideSet(peptides, pcssRunner.getNegativeKeyword())
        self.totalPositiveCount = len(self.positivePeptides)
        self.totalNegativeCount = len(self.negativePeptides)

        fraction = pcssRunner.pcssConfig["jackknife_fraction"]
        self.testSetPositiveCount = math.floor(float(self.totalPositiveCount) * float(fraction))
        self.trainingSetPositiveCount = self.totalPositiveCount - self.testSetPositiveCount
        self.trainingSetNegativeCount = self.trainingSetPositiveCount
        self.testSetNegativeCount = self.totalNegativeCount - self.trainingSetNegativeCount
        self.pcssRunner = pcssRunner
        self.validateCounts()

        print "got %s positive, %s negative,  %s training positive, %s test positive, %s training negative, %s test negatives" % (self.totalPositiveCount, self.totalNegativeCount, 
                                                                                                                                  self.trainingSetPositiveCount, 
                                                                                                                                  self.testSetPositiveCount, 
                                                                                                                                  self.trainingSetNegativeCount, 
                                                                                                                                  self.testSetNegativeCount)
    def validateCounts(self):
        if (self.testSetPositiveCount < 1):
            raise pcssErrors.PcssGlobalException("Positive test set count is %s (should be greater than 0)")
        if (self.testSetPositiveCount > self.testSetNegativeCount):
            raise pcssErrors.PcssGlobalException("Test set should have more negatives than positives")

        if (self.testSetPositiveCount + self.trainingSetPositiveCount != self.totalPositiveCount):
            raise pcssErrors.PcssGlobalException("Positive Training Set (%s) and Positive Test Set (%s) do not add up to total positives (%s)"
                                                 % (self.trainingSetPositiveCount, self.testSetPositiveCount, self.totalPositiveCount))
        if (self.testSetNegativeCount + self.trainingSetNegativeCount != self.totalNegativeCount):
            raise pcssErrors.PcssGlobalException("Negative Training Set (%s) and Negative Test Set (%s) do not add up to total negatives (%s)"
                                                 % (self.trainingSetNegativeCount, self.testSetNegativeCount, self.totalNegativeCount))
        

    def makeTrainingAndTestSets(self):
        self.makePositiveTestSet()
        self.makeNegativeTestSet()
        self.makePositiveTrainingSet()
        self.makeNegativeTrainingSet()
        
    def makePeptideSet(self, peptides, statusType):
        peptideSet = []
        for peptide in peptides:
            if (peptide.getAttributeOutputString("status") == statusType):
                peptideSet.append(peptide)
        return peptideSet

    def getSample(self, peptides, count):
        if (count > len(peptides)):
            raise pcssErrors.PcssGlobalException("getSample(): tried to sample %s peptides but there are only %s peptides in the pool" 
                                                 % (count, len(peptides)))
        makeRandomSample = self.pcssRunner.internalConfig["make_random_test_set"]
        if (makeRandomSample):
            print "RANDOM SAMPLE"
            return random.sample(peptides, int(count))
        else:
            print "NON RANDOM SAMPLE" # -- make internal config interpolation and test
            return peptides[0:int(count)]

    def makePositiveTestSet(self):
        self.positiveTestSet = self.getSample(self.positivePeptides, self.testSetPositiveCount)
        
    def makeNegativeTestSet(self):
        self.negativeTestSet = self.getSample(self.negativePeptides, self.testSetNegativeCount)

    def makePositiveTrainingSet(self):
        self.positiveTrainingSet = self.getRemainingPeptides(self.positivePeptides, self.positiveTestSet)

    def makeNegativeTrainingSet(self):
        self.negativeTrainingSet = self.getRemainingPeptides(self.negativePeptides, self.negativeTestSet)

    def getRemainingPeptides(self, allPeptides, peptidesSoFar):
        remainingPeptides = []
        for peptide in allPeptides:
            if peptide not in peptidesSoFar:
                remainingPeptides.append(peptide)
        return remainingPeptides
        

class ClassifySvm:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.PeptideScoreTuple = collections.namedtuple('peptideScore', ['peptide', 'score'])

        self.pstList = []

    def setProteins(self, proteins):
        self.proteins = proteins
        peptides = []
        for protein in proteins:
            if (not protein.hasErrors()):
                nextPeptides = protein.peptides.values()
                for peptide in nextPeptides:
                    peptides.append(peptide)
        self.setPeptides(peptides)

    def setPeptides(self, peptides):
        self.peptides = peptides
        
    def writeClassificationFile(self):
        classificationFileName = self.getSvmInputFile()
        classificationFh = open(classificationFileName, 'w')
        i = 1
        for peptide in self.peptides:
            print "write peptide %s count %s" % (peptide.startPosition, i)
            i += 1
            nextLine = peptide.makeSvmFileLine()
            classificationFh.write("%s %s\n" % ("0", nextLine))
        classificationFh.close()

    def classifySvm(self):
        svmCommandName = self.pcssRunner.internalConfig['svm_classify_command']
        classificationFileName = self.getSvmInputFile()

        if (not os.path.exists(classificationFileName)):
            raise pcssErrors.PcssGlobalException("Did not find test set file in expected location -- searched for\n%s" % classificationFileName)

        scoreFileName = self.getClassifyOutputFile()

        modelFile = self.getSvmModelFile()

        svmOutput = self.pcssRunner.pdh.runSubprocess([svmCommandName, classificationFileName, modelFile, scoreFileName])
        
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
            pst = self.PeptideScoreTuple(peptide, score)
            self.pstList.append(pst)

class BenchmarkPeptideScoreTupleList:
    def __init__(self, pcssRunner):
        self.pstList = []
        self.pcssRunner = pcssRunner
        self.BenchmarkPeptideScoreTuple = collections.namedtuple('peptideScore', ['peptide', 'score', 'positiveCount', 'negativeCount', 'tpr', 'fpr'])

    def addPst(self, peptideScoreTuple):

        self.pstList.append(peptideScoreTuple)

    def getPeptideStatusCount(self, pstList, status):
        totalCount = 0
        for i in pstList:
            if (i.peptide.getAttributeOutputString("status") == status):
                totalCount += 1
        return totalCount

    def finalize(self):
        sortedPstList = sorted(self.pstList, key=lambda pst: pst.score)
        self.benchmarkPstList = []
        print "have %s peptides" % len(sortedPstList)
        self.totalPositiveCount = self.getPeptideStatusCount(sortedPstList, self.pcssRunner.getPositiveKeyword())
    
        self.totalNegativeCount = self.getPeptideStatusCount(sortedPstList, self.pcssRunner.getNegativeKeyword())
        print "%s positive peptides and %s negative peptides" % (self.totalPositiveCount, self.totalNegativeCount)
        currentPositiveCount = 0
        currentNegativeCount = 0

        for pst in sortedPstList:
            if (pst.peptide.getAttributeOutputString("status") == self.pcssRunner.getPositiveKeyword()):
                currentPositiveCount += 1
            elif (pst.peptide.getAttributeOutputString("status") == self.pcssRunner.getNegativeKeyword()):
                currentNegativeCount += 1
                
            bpst = self.BenchmarkPeptideScoreTuple(pst.peptide, pst.score, currentPositiveCount, currentNegativeCount, 
                                                   self.getRate(currentPositiveCount, self.totalPositiveCount), 
                                                   self.getRate(currentNegativeCount, self.totalNegativeCount))
            self.benchmarkPstList.append(bpst)

    def getRate(self, currentCount, totalCount):
        return float(currentCount) / float(totalCount)

    def getBenchmarkTuple(self, i):
        return self.benchmarkPstList[i]

class ApplicationSvm(ClassifySvm):

    def addScoresToPeptides(self):
        benchmarkResultsFile =  self.pcssRunner.pcssConfig["svm_benchmark_file"]

        self.br = BenchmarkResults()
        self.br.readBenchmarkFile(benchmarkResultsFile)
        print "adding scores to peptides"
        for pst in self.pstList:
            st = self.br.getClosestScoreTuple(pst.score)
            pst.peptide.addStringAttribute('svm_score', st.score)
            pst.peptide.addStringAttribute('svm_fpr', round(st.fpr, 3))
            pst.peptide.addStringAttribute('svm_tpr', round(st.tpr, 3))

    def getSvmInputFile(self):
        return self.pcssRunner.pdh.getSvmApplicationSetFile()
    
    def getClassifyOutputFile(self):
        return self.pcssRunner.pdh.getSvmApplicationOutputFile()

    def getSvmModelFile(self):
        return self.pcssRunner.pcssConfig["svm_model_file"]


class TestSvm(ClassifySvm):
    def getSvmInputFile(self):
        return self.pcssRunner.pdh.getSvmTestSetFile()
    
    def getClassifyOutputFile(self):
        return self.pcssRunner.pdh.getSvmTestOutputFile()

    def getSvmModelFile(self):
        modelFileName = self.pcssRunner.pdh.getSvmNewModelFile()
        if (not os.path.exists(modelFileName)):
            raise pcssErrors.PcssGlobalException("Could not find newly created model file for test svm (looked for %s)" % modelFileName)
        return modelFileName
                                                

    def getBenchmarkScoreTupleList(self):
        bpstList = BenchmarkPeptideScoreTupleList(self.pcssRunner)
        for pst in self.pstList:
            bpstList.addPst(pst)
        bpstList.finalize()
        return bpstList

class SvmFeatureHandler:
    def __init__(self):
        self.featureNumber = 1

    def processEmptyFeature(self, peptide, feature):
        
        self.featureNumber += feature.getEmptyFeatureOffset(peptide)
        
    def getFeatureNumber(self):
        return self.featureNumber

    def processFeature(self, increment):
        self.featureNumber += increment
