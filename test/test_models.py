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
import pcssIO
import pcssTests
from Bio import PDB

class TestModels(pcssTests.PcssTest):

    def processModelException(self, exceptionCode, function):

        peptide = self.proteins[0].peptides.values()[0]
        self.processFeatureException(peptide, "dssp_structure", exceptionCode, function)

    def setupSpecificTest(self):
        self.runner = pcssTools.AnnotationRunner(self.pcssConfig)
        spi = pcssIO.ScanPeptideImporter(self.runner)
        self.proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        currentModelFile = self.runner.pdh.getFullModelFileFromId("741bc8ce184702f143409644b7a6f690")

    def addModelsToTestProtein(self):
        pcssProtein = self.createTestProtein()
        pcssProtein.addModels(self.modelTable)
        return pcssProtein

    def createTestProtein(self):
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)
        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ", self.proteins)
        return pcssProtein

    def test_missing_source_file(self):
        pcssProtein = self.addModelsToTestProtein()
        bestModel = pcssProtein.getRankedModels()[0]
        bestModel.setAttribute("model_id", "fake")
        self.processModelException("peptide_error_no_source_model", pcssProtein.processDssp)

    def test_invalid_model_range(self):
        self.pcssConfig['model_table_column_file'] = "testInput/dsspErrors/modelColumnOrderBadRange.txt"            
        pcssProtein = self.createTestProtein()
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            pcssProtein.addModels(self.modelTable)
        self.handleTestException(pge)

    def test_dssp_error(self):
        pcssProtein = self.addModelsToTestProtein()
        self.pcssConfig['dssp_executable'] = "fake"
        self.processModelException("peptide_error_dssp_error", pcssProtein.processDssp)

    def test_dssp_peptide_mismatch(self):
        pcssProtein = self.addModelsToTestProtein()
        pcssProtein.peptides[2].sequence = "FAKEFAKE"
        self.processModelException("peptide_error_dssp_mismatch", pcssProtein.processDssp)

    def test_column_count_mismatch(self):
        self.pcssConfig['model_table_column_file'] = "testInput/dsspErrors/modelColumnOrderShort.txt"    
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            pcssModels.PcssModelTable(self.runner, modelColumns)
        self.handleTestException(pge)

    def test_read_model_table(self):

        pcssProtein = self.addModelsToTestProtein()

        modelTableSequence = self.modelTable.getPcssModelSequence(pcssProtein.modbaseSequenceId)            
        models = modelTableSequence.getModels()
        self.assertEquals(len(models), 4)

        pcssModel = modelTableSequence.getModel("741bc8ce184702f143409644b7a6f690")
        self.assertTrue(pcssModel.getId(),"741bc8ce184702f143409644b7a6f690") 

        self.assertEquals(float(pcssModel.getAttributeValue("no35")), 0.991304)

        pcssProtein.processDssp()
        rankedModels = pcssProtein.getRankedModels()
        self.assertEquals(rankedModels[0].getId(), pcssModel.getId())

        peptide = pcssProtein.peptides[2]
        self.assertEquals(peptide.bestModel.getId(), pcssModel.getId())

        with self.assertRaises(pcssErrors.PcssGlobalException) as pge:
            pcssModel.getAttributeValue("fake")
        self.handleTestException(pge)

        self.assertTrue(os.path.exists(self.runner.pdh.getFullModelFile(pcssModel)))

        self.assertEquals(self.runner.getModelUrl("MYMODELID"), "http://salilab.org/modbase/search?modelID=MYMODELID&displaymode=moddetail")

        self.assertEquals(peptide.attributes["dssp_structure"].getValueString(), "AAAAAAAA")
        self.assertEquals(peptide.attributes["dssp_accessibility"].getValueString(), 
                          "0.589, 0.234, 0.763, 0.515, 0.189, 0.162, 0.658, 0.374")

if __name__ == '__main__':
    unittest.main()


