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
            

class CompleteSvmGenerator:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner

    def createSvmModel(self, peptides):
        self.trainingSvm = TrainingSvm(self.pcssRunner)
        self.trainingSvm.setPeptides(peptides)
        modelFile = self.pcssRunner.pdh.getUserModelFileName()
        self.trainingSvm.writePeptidesToFile(modelFile)
        

class LeaveOneOutBenchmarker:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.trainingSvm = TrainingSvm(pcssRunner)
        self.testSvm = TestSvm(pcssRunner)  
        self.currentPeptidePosition = 0
        self.looTsr = TestSetResult(pcssRunner)

    def createTrainingAndTestSets(self, peptides):
        trainingPeptideList = []
        if (self.currentPeptidePosition >= len(peptides)):
            msg = "Error: Leave one out benchmarker internal peptide counter (%s) must be smaller than input peptide set count (%s)" % (self.currentPeptidePosition,
                                                                                                                                        len(peptides))
            raise pcssErrors.PcssGlobalException(msg)
        for (i, peptide) in enumerate(peptides):
            if (i == self.currentPeptidePosition):
                self.testSvm.setPeptides([peptide])
                print "next test set peptide position %s status %s" % (peptide.startPosition, peptide.getAttributeOutputString("status"))
            else:
                trainingPeptideList.append(peptide)
        self.trainingSvm.setPeptides(trainingPeptideList)
        self.trainingSvm.writeTrainingSetFile()
        self.testSvm.writeClassificationFile()
        self.currentPeptidePosition += 1
        
    def trainAndApplyModel(self):
        
        self.trainingSvm.trainModel()

        if (len(self.testSvm.peptides) > 1):
            raise pcssErrors.PcssGlobalException("Error: Leave One Out Benchmarker has test set greater than size 1 (%s total)" % len(self.testSvm.peptides))
        
        self.testSvm.classifySvm()

    def readBenchmarkResults(self):
        self.testSvm.readResultFile()
        
        pstList = self.testSvm.getPstList()
        
        self.looTsr.addPst(pstList[0])
        
    def processAllResults(self):
        self.looTsr.finalize()

        resultFile = self.pcssRunner.pdh.getLeaveOneOutResultFileName()
        resultFh = open(resultFile, 'w')
        size = self.looTsr.getSize()
        for i in range(size):
            nextBpst = self.looTsr.getBenchmarkTuple(i)
            #add petpides start position and modbase seq id
            outputList = [nextBpst.fpr, nextBpst.tpr, nextBpst.score,  nextBpst.peptide.getAttributeOutputString("status"), 
                          nextBpst.negativeCount, nextBpst.positiveCount]
            resultFh.write("%s\n" % "\t".join(str(x) for x in outputList))

class SvmBenchmarker:
    def __init__(self, pcssRunner):
        self.pcssRunner = pcssRunner
        self.trainingSvm = TrainingSvm(pcssRunner)
        self.testSvm = TestSvm(pcssRunner)  
        self.testSetResultTracker = TestSetResultTracker(self.pcssRunner)
        
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
        
        testSetResult = self.testSvm.getTestSetResult()
        #might want to write to file if not doing everything in memory
        self.testSetResultTracker.addTestSetResult(testSetResult)

    def processAllResults(self):

        self.testSetResultTracker.finalize()
        for i in range(self.testSetResultTracker.getTprCount()):
            nextTuple = self.testSetResultTracker.getBenchmarkTuple(i)
            print "next run: tpr %s fpr %s score  %s stddev %s" % (nextTuple.tpr, nextTuple.fpr, nextTuple.score, nextTuple.fprStdev)
        
        self.testSetResultTracker.writeResultFile(self.pcssRunner.pdh.getFullBenchmarkResultFileName())
        
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
    
    def writeFullPeptideModelFile(self):
        fullPeptideModelFileName = self.runner.pdh.getUserCreatedModelFileName()
        self.writePeptidesToFile(fullPeptideModelFileName)

    def writeTrainingSetFile(self):
        trainingSetFileName  = self.runner.pdh.getSvmTrainingSetFile()
        self.writePeptidesToFile(trainingSetFileName)
    
    def writePeptidesToFile(self, fileName):
        
        trainingSetFh = open(fileName, 'w')
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
        self.peptides = []        
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
        print "set peptides; have %s total" % len(self.peptides)
        
    def writeClassificationFile(self):
        classificationFileName = self.getSvmInputFile()
        classificationFh = open(classificationFileName, 'w')
        i = 1
        for peptide in self.peptides:
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
        self.pstList = []
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

    def getPstList(self):
        return self.pstList

class TestSetResult:
    def __init__(self, pcssRunner):
        self.pstList = []
        self.pcssRunner = pcssRunner
        self.setupTuple()

    def setupTuple(self):
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

    def getFinalPositiveCount(self):
        return self.benchmarkPstList[-1].positiveCount

    def getFinalNegativeCount(self):
        return self.benchmarkPstList[-1].negativeCount

    def getIncrementedTprTuples(self):
        currentTpr = 0.0
        tuplePositions = []
        for bpst in self.benchmarkPstList:
            if (bpst.tpr > currentTpr):
                currentTpr = bpst.tpr
                tuplePositions.append(bpst)
        return tuplePositions

    def getRate(self, currentCount, totalCount):
        return float(currentCount) / float(totalCount)

    def getBenchmarkTuple(self, i):
        return self.benchmarkPstList[i]

    def getSize(self):
        return len(self.benchmarkPstList)

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

class TestSetResultTracker:

    def __init__(self, runner):
        self.pcssRunner = runner
        self.TestSetTprAveragesTuple = collections.namedtuple('singleTpr', ['tpr', 'fpr', 'score', 'fprStdev'])
        self.allTestSetResults = []
        
    def addTestSetResult(self, pstList):
        
        self.allTestSetResults.append(pstList)

    def validateTestSetResult(self, tsr):
        if (tsr.getFinalPositiveCount() != self.referencePositiveCount):
            raise pcssErrors.PcssGlobalException("Error: test set did not have same number of positives (%s) as the reference (%s)" % (tsr.getFinalPositiveCount(),
                                                                                                                                       self.referencePositiveCount))
        if (tsr.getFinalNegativeCount() != self.referenceNegativeCount):
            raise pcssErrors.PcssGlobalException("Error: test set did not have same number of negatives (%s) as the reference (%s)" % (tsr.getFinalNegativeCount(),
                                                                                                                                       self.referenceNegativeCount))
        
    def finalize(self):
        tprCountsToTuples = {}
        self.referencePositiveCount = self.allTestSetResults[0].getFinalPositiveCount()
        self.referenceNegativeCount = self.allTestSetResults[0].getFinalNegativeCount()
        for testSetResult in self.allTestSetResults:
            self.validateTestSetResult(testSetResult)
            tprTuples = testSetResult.getIncrementedTprTuples()
            for tprTuple in tprTuples:
                if (not tprTuple.tpr in tprCountsToTuples):
                    tprCountsToTuples[tprTuple.tpr] = []
                tprCountsToTuples[tprTuple.tpr].append(tprTuple)
            #assert same FPR?
        self.testSetTprAveragesList = []
        for tpr, bpstList in sorted(tprCountsToTuples.iteritems()):

            fprAverage = self.average(list(x.fpr for x in bpstList))
            fprStdev = self.stddev(list(x.fpr for x in bpstList))
            scoreAverage = self.average(list(x.score for x in bpstList))
            self.testSetTprAveragesList.append(self.TestSetTprAveragesTuple(tpr, fprAverage, scoreAverage, fprStdev))

    def getBenchmarkTuple(self, i):
        return self.testSetTprAveragesList[i]

    def getRunCount(self):
        return len(self.allTestSetResults)

    def getTprCount(self):
        return len(self.testSetTprAveragesList)

    def average(self, list):
        average = float(sum(list)) / float(len(list))
        return average

    def stddev(self, list):
        stddev = math.sqrt((float(sum(x * x for x in list)) / float(len(list))) - (float(self.average(list)) * float(self.average(list))))
        return stddev
    

    def writeResultFile(self, resultFileName):
        resultFh = open(resultFileName, 'w')
        initialOutputList = ["0", "0", "N/A", "N/A"]
        finalOutputList = ["1", "1", "N/A", "N/A"]
        resultFh.write("%s\n" % "\t".join(initialOutputList))
        for averagesTuple in self.testSetTprAveragesList:
            outputList = [averagesTuple.fpr, averagesTuple.tpr, averagesTuple.score, averagesTuple.fprStdev]
            resultFh.write("%s\n" % "\t".join(str(round(x, 3)) for x in outputList))
        resultFh.write("%s\n" % "\t".join(finalOutputList))


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
                                                
    def getTestSetResult(self):
        tsr = TestSetResult(self.pcssRunner)
        for pst in self.pstList:
            tsr.addPst(pst)
        tsr.finalize()
        return tsr
    

class SvmFeatureHandler:
    def __init__(self):
        self.featureNumber = 1

    def processEmptyFeature(self, peptideLength, feature):
        
        self.featureNumber += feature.getEmptyFeatureOffset(peptideLength)
        
    def getFeatureNumber(self):
        return self.featureNumber

    def processFeature(self, increment):
        self.featureNumber += increment

    def finalizeFeature(self, peptide, feature, referencePeptideLength):
        if (peptide.getPeptideLength() > referencePeptideLength):
            raise pcssErrors.PcssGlobalException("Peptide %s has length of %s which is greater than reference %s" % (peptide.startPosition, 
                                                                                                                     peptide.getPeptideLength(),
                                                                                                                     referencePeptideLength()))
        lengthDifference = referencePeptideLength - peptide.getPeptideLength() 
        multiplier = feature.getFeatureLength()
        self.featureNumber += (lengthDifference * multiplier)
