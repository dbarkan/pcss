import sys
import os
import re
import subprocess
import configobj
import pcssTools
import pcssPeptide
import pcssErrors
import pcssIO
import pcssModels
import pcssFeatures
import pcssFeatureHandlers
import logging

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)

class TempPcssRunner:

    def __init__(self, configFile):

        self.pcssConfig = configobj.ConfigObj(configFile)
        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.pcssConfig['fasta_file'])
        
    def runStructureFeatures(self):
        
        try:
            modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
            modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)
            for protein in self.proteins:
                protein.addModels(modelTable)
                protein.processDssp()
        except pcssErrors.PcssShutilError as e:
            print str(e.function)
            print e.args
            print e.ioError
        except pcssErrors.PcssGlobalException as e:
            print e.msg

    def runSequenceFeatures(self):
        
        try:
            
            disopredFileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.runner.pdh)
            disopredReader = pcssFeatureHandlers.DisopredReader(disopredFileHandler)
            disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(disopredFileHandler)

            for protein in self.proteins:
                protein.processDisopred(disopredReader, disopredRunner)


            psipredFileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.runner.pdh)
            psipredReader = pcssFeatureHandlers.PsipredReader(psipredFileHandler)
            psipredRunner = pcssFeatureHandlers.SequenceFeatureRunner(psipredFileHandler)

            for protein in self.proteins:
                protein.processPsipred(psipredReader, psipredRunner)
                
            print self.proteins[0].getOutputString()

        except pcssErrors.PcssGlobalException as e:
            print "got Pcss Global error"
            print "Message: %s %s" % (e.msg, e.expr)

    def writeOutput(self):
        afw = pcssIO.AnnotationFileWriter(self.runner)
        afw.writeAllOutput(self.proteins)
        
if __name__ == '__main__':

    configFile = sys.argv[1]
    tempRunner = TempPcssRunner(configFile)
    tempRunner.runSequenceFeatures()
    tempRunner.runStructureFeatures()
    tempRunner.writeOutput()
