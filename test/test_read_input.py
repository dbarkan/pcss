import unittest
import configobj
import pcssTools
import pcssPeptide
import pcssModels
import logging
import pcssErrors
import pcssIO
import pcssFeatureHandlers
import pcssFeatures
import os
import sys

logging.basicConfig(stream=sys.stdout)
logging.root.setLevel(logging.DEBUG)


class TestReadInput(unittest.TestCase):

    def getProtein(self, proteinId):
        for protein in self.proteins:
            if (protein.modbaseSequenceId == proteinId):
                return protein

    def setUp(self):
        configFile = "testConfig/testPcssConfig.txt"
        configSpecFile = "testConfig/testConfigSpec.txt"
        self.pcssConfig = configobj.ConfigObj(configFile, configspec=configSpecFile)

        self.runner = pcssTools.PcssRunner(self.pcssConfig)
        self.spi = pcssIO.ScanPeptideImporter(self.runner)

          
    def test_bad_annotation_file(self):
        annotationFileName = os.path.join(self.pcssConfig["home_test_directory"],
                                          "testInput/ioErrors/missingColumnsFile.txt")
        reader = pcssIO.AnnotationFileReader(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readAnnotationFile(annotationFileName)
        print e.exception.msg

        annotationFileName = os.path.join(self.pcssConfig["home_test_directory"], "fake")
        self.assertRaises(pcssErrors.PcssGlobalException, reader.readAnnotationFile, annotationFileName)
        
        annotationFileName = os.path.join(self.pcssConfig["home_test_directory"],
                                          "testInput/ioErrors/missingInputColumn.txt")
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readAnnotationFile(annotationFileName)
        print e.exception.msg

        annotationFileName = os.path.join(self.pcssConfig["home_test_directory"],
                                          "testInput/ioErrors/extraInputColumn.txt")
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readAnnotationFile(annotationFileName)
        print e.exception.msg


    def test_empty_rules_file(self):
        self.pcssConfig["rules_file"] = "testInput/emptyPeptideRulesFile"

        spi = pcssIO.ScanPeptideImporter(self.runner)
        proteins = spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        firstProtein = proteins[0]
        self.assertEquals(248, len(firstProtein.peptides.values()))

    def test_bad_protein_attribute(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        badAtt = self.runner.pfa.getColumnSortedInputAttributes()[0]
        badAtt.name = "seq_id_fake"
        
        afw = pcssIO.AnnotationFileWriter(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            afw.writeAllOutput(self.proteins)
        print e.exception.msg
    
    def test_bad_peptide_attribute(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        startAtt = self.runner.pfa.getAttribute('peptide_start')
        startAtt.name = "peptide_start_fake"
        afw = pcssIO.AnnotationFileWriter(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            afw.writeAllOutput(self.proteins)
        print e.exception.msg
    
    def test_bad_get_protein_attribute(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            self.proteins[0].setStringAttribute("fake", "fakeValue")
        print e.exception.msg
    
    def test_read_attributes(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEquals(self.runner.pfa.getColumnSortedInputAttributes()[-1].name, "peptide_errors")

    def test_write_normal_output(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)
        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ")

        disopredFileHandler = pcssFeatureHandlers.DisopredFileHandler(self.pcssConfig, self.runner.pdh)
        disopredReader = pcssFeatureHandlers.DisopredReader(disopredFileHandler)
        disopredRunner = pcssFeatureHandlers.SequenceFeatureRunner(disopredFileHandler)
        pcssProtein.processDisopred(disopredReader, disopredRunner)

        psipredFileHandler = pcssFeatureHandlers.PsipredFileHandler(self.pcssConfig, self.runner.pdh)
        psipredReader = pcssFeatureHandlers.PsipredReader(psipredFileHandler)
        psipredRunner = pcssFeatureHandlers.SequenceFeatureRunner(psipredFileHandler)
        pcssProtein.processPsipred(psipredReader, psipredRunner)

        pcssProtein.addModels(self.modelTable)        
        pcssProtein.processDssp()
        rankedModels = pcssProtein.getRankedModels()
        try:
            afw = pcssIO.AnnotationFileWriter(self.runner)
            afw.writeAllOutput(self.proteins)
        except Exception as e:
            print e
            print e.msg

        annotationFileName = self.runner.pdh.getFullOutputFile("annotationOutput.txt")
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile(annotationFileName)
        newProteins = reader.getProteins()


        for newProtein in newProteins:

            oldProtein = self.getProtein(newProtein.modbaseSequenceId)
            self.assertTrue(oldProtein.isEqual(newProtein))
        

    def test_read_feature_error(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/ioErrors/annotationOutputFeatureError.txt")
        firstProtein = reader.getProteins()[0]
        firstPeptide = firstProtein.peptides.values()[0]
        self.assertEqual("peptide_error_no_source_model", firstPeptide.getAttributeOutputString("dssp_structure"))

    def test_no_peptides_parsed(self):
        self.runner.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["home_test_directory"],
                                                            "testInput/ioErrors/noPeptidesParsedFasta.txt")

        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEqual(self.proteins[0].getAttributeOutputString("protein_errors"), "no_peptides_parsed")
        afw = pcssIO.AnnotationFileWriter(self.runner)
        self.runner.pfa.setAllOptional()

    def test_annotation_no_sequences(self):
        self.runner.pcssConfig['fasta_file'] = os.path.join(self.pcssConfig["home_test_directory"],
                                                            "testInput/ioErrors/annotationMissingSequence.txt")
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/ioErrors/normalAnnotationFile.txt")
        self.assertRaises(pcssErrors.PcssGlobalException, reader.readProteinSequences, self.runner.pcssConfig['fasta_file'])

    def test_scan_peptides(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEqual(self.proteins[0].modbaseSequenceId, "76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        self.assertEqual(self.proteins[0].uniprotId, "P62258")
        self.assertEqual(self.proteins[0].peptides.values()[0].startPosition, 2)
        self.assertEqual(self.proteins[0].peptides.values()[0].endPosition, 9)
        self.assertEqual(self.proteins[0].peptides.values()[0].sequence, "DREDLVYQ")
        self.assertEqual(len(self.proteins[0].peptides.values()), 19)
        

    def test_defined_peptide_mismatch(self):
        self.runner.pcssConfig['fasta_file'] = "testInput/ioErrors/peptideMismatchFasta.txt"
        dpi = pcssIO.DefinedPeptideImporter(self.runner)
        self.assertRaises(pcssErrors.PcssGlobalException, dpi.readInputFile, self.runner.pcssConfig['fasta_file'])
        
    def test_defined_peptides(self):
        self.runner.pcssConfig['fasta_file'] = "testInput/inputSequenceDefined.txt"
        dpi = pcssIO.DefinedPeptideImporter(self.runner)
    
        self.proteins = dpi.readInputFile(self.runner.pcssConfig['fasta_file'])
        
        self.assertEqual(self.proteins[0].modbaseSequenceId, "76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        self.assertEqual(self.proteins[0].uniprotId, "P62258")
        self.assertEqual(self.proteins[0].peptides.values()[0].startPosition, 2)
        self.assertEqual(self.proteins[0].peptides.values()[0].endPosition, 9)
        self.assertEqual(self.proteins[0].peptides.values()[0].sequence, "DREDLVYQ")

    def test_bad_status_fasta(self):
        self.runner.pcssConfig['fasta_file'] = "testInput/ioErrors/badStatusFasta.txt"
        dpi = pcssIO.DefinedPeptideImporter(self.runner)
        
        self.assertRaises(pcssErrors.PcssGlobalException, dpi.readInputFile, self.runner.pcssConfig['fasta_file'])

        self.runner.pcssConfig['fasta_file'] = "testInput/ioErrors/noPeptidesFasta.txt"
        self.assertRaises(pcssErrors.PcssGlobalException, dpi.readInputFile, self.runner.pcssConfig['fasta_file'])
    
        self.runner.pcssConfig['fasta_file'] = "testInput/ioErrors/onlyMismatchFasta.txt"
        self.assertRaises(pcssErrors.PcssGlobalException, dpi.readInputFile, self.runner.pcssConfig['fasta_file'])

    def test_no_annotation_proteins(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        self.assertRaises(pcssErrors.PcssGlobalException, reader.readAnnotationFile, "testInput/ioErrors/noAnnotationProteins.txt")
        
    def read_write_annotation_file(self, inputFile):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile(inputFile)
        proteins = reader.getProteins()
        afw = pcssIO.AnnotationFileWriter(self.runner)
        afw.writeAllOutput(proteins)
        self.compareFiles(inputFile, 
                          self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]), True)
        
    def test_read_write_normal_annotation_file(self):    
        self.read_write_annotation_file("testInput/ioErrors/normalAnnotationFile.txt")

    def test_read_write_feature_error_file(self):
        self.read_write_annotation_file("testInput/ioErrors/annotationOutputFeatureError.txt")
        

    def compareFiles(self, firstFile, secondFile, sortLines=False):
        firstReader = pcssTools.PcssFileReader(firstFile)
        secondReader = pcssTools.PcssFileReader(secondFile)
        if (sortLines):
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())
        else:
            firstLines =  sorted(firstReader.getLines())
            secondLines = sorted(secondReader.getLines())

        for (i, firstLine) in enumerate(firstLines):

            secondLine = secondLines[i]
            self.assertEquals(firstLine, secondLine)


if __name__ == '__main__':
    unittest.main()
