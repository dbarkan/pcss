import unittest
import os
import configobj
import pcssTools
import pcssPeptide
import pcssModels
import pcssFeatures
import pcssFeatureHandlers
import pcssErrors
import runpy
from Bio import PDB

class TestModels(unittest.TestCase):
    
    def setupPcssModelTest(self):
        configFile = "testConfig/testPcssConfig.txt"

        self.pcssConfig = configobj.ConfigObj(configFile)

        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        
        spi = pcssPeptide.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        currentModelFile = self.runner.pdh.getFullModelFileFromId("741bc8ce184702f143409644b7a6f690")
        if (os.path.exists(currentModelFile)):
            os.remove(currentModelFile)


    def getProtein(self, proteinId):
        for protein in self.proteins:

            if (protein.modbaseSequenceId == proteinId):
                return protein

    def test_missing_source_file(self):

        self.setupPcssModelTest()
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)

        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        pcssProtein.addModels(self.modelTable)
        bestModel = pcssProtein.getRankedModels()[0]
        print "best model id: %s" % bestModel.getId()
        bestModel.setAttribute("model_id", "fake")

        with self.assertRaises(pcssErrors.FeatureException) as fe:
            pcssProtein.processDssp(self.modelTable)
        self.assertTrue(fe.exception.pcssPeptide, not None)                
        self.assertEquals(fe.exception.pcssProtein.modbaseSequenceId, self.proteins[0].modbaseSequenceId)
        
    def test_column_count_mismatch(self):
        self.setupPcssModelTest()
        self.pcssConfig['model_table_column_file'] = "testInput/dsspErrors/modelColumnOrderShort.txt"    
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.assertRaises(pcssErrors.PcssGlobalException, pcssModels.PcssModelTable, self.runner, modelColumns)

    def test_invalid_model_range(self):
        self.setupPcssModelTest()
        self.pcssConfig['model_table_column_file'] = "testInput/dsspErrors/modelColumnOrderBadRange.txt"    
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)

        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        self.assertRaises(pcssErrors.PcssGlobalException, pcssProtein.addModels, self.modelTable)


    def test_dssp_error(self):
        self.setupPcssModelTest()
        self.pcssConfig['dssp_executable'] = "fake"
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)

        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        pcssProtein.addModels(self.modelTable)
        with self.assertRaises(pcssErrors.FeatureException) as fe:
            pcssProtein.processDssp(self.modelTable)
        self.assertTrue(fe.exception.pcssPeptide, not None)                
        self.assertEquals(fe.exception.pcssProtein.modbaseSequenceId, self.proteins[0].modbaseSequenceId)
        print fe.exception.msg

    def test_dssp_peptide_mismatch(self):
        self.setupPcssModelTest()
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)

        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        pcssProtein.addModels(self.modelTable)
        pcssProtein.peptides[2].sequence = "FAKEFAKE"
        with self.assertRaises(pcssErrors.FeatureException) as fe:
            pcssProtein.processDssp(self.modelTable)
        self.assertTrue(fe.exception.pcssPeptide, not None)                
        self.assertEquals(fe.exception.pcssProtein.modbaseSequenceId, self.proteins[0].modbaseSequenceId)
        print fe.exception.msg

    def test_read_model_table(self):
        try:
            self.setupPcssModelTest()
            modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
            self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)
            
            pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")
            modelTableSequence = self.modelTable.getPcssModelSequence(pcssProtein.modbaseSequenceId)            

            models = modelTableSequence.getModels()
            self.assertEquals(len(models), 4)

            pcssModel = modelTableSequence.getModel("741bc8ce184702f143409644b7a6f690")
            self.assertTrue(pcssModel, not None)

            self.assertEquals(float(pcssModel.getAttributeValue("no35")), 0.991304)
            pcssProtein.addModels(self.modelTable)        
            pcssProtein.processDssp(self.modelTable)
            rankedModels = pcssProtein.getRankedModels()
            self.assertEquals(rankedModels[0].getId(), pcssModel.getId())

            peptide = pcssProtein.peptides[2]
            self.assertEquals(peptide.bestModel.getId(), pcssModel.getId())
            
            self.assertRaises(pcssErrors.PcssGlobalException, pcssModel.getAttributeValue, "fake")
            self.assertTrue(os.path.exists(self.runner.pdh.getFullModelFile(pcssModel)))

            self.assertEquals(peptide.features["dssp_structure_feature"].getValueString(), "AAAAAAAA")
            self.assertEquals(peptide.features["dssp_accessibility_feature"].getValueString(), 
                              "0.589, 0.234, 0.763, 0.515, 0.189, 0.162, 0.658, 0.374")

        except pcssErrors.PcssGlobalException as e:
            print e.msg

        except pcssErrors.FeatureException as e:
            print e.msg
            print e.pcssPeptide.sequence
            print e.pcssPeptide.startPosition


if __name__ == '__main__':
    unittest.main()


