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
import pcssTests

class TestReadInput(pcssTests.PcssTest):
    def setupSpecificTest(self):
        
        self.runner = pcssTools.ModelRunner(self.pcssConfig)
        self.spi = pcssIO.ScanPeptideImporter(self.runner)
        self.errorDirName = "ioErrors"
        self.testName = "io"
          
    def processBadAnnotationFileTest(self, annotationFileName):
        reader = pcssIO.AnnotationFileReader(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readAnnotationFile(annotationFileName)
        self.handleTestException(e)
        print "annotation file name: %s exception msg: %s" % (annotationFileName, e.exception.msg)

    def test_bad_annotation_file(self):
        self.processBadAnnotationFileTest(self.getErrorInputFile("missingColumnsFile.txt"))
        self.processBadAnnotationFileTest("fake")
        self.processBadAnnotationFileTest(self.getErrorInputFile("missingInputColumn.txt"))
        self.processBadAnnotationFileTest(self.getErrorInputFile("extraInputColumn.txt"))
        
    def test_empty_rules_file(self):
        self.pcssConfig["rules_file"] = self.getErrorInputFile("emptyPeptideRulesFile")
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
        self.handleTestException(e)
    
    def test_bad_peptide_attribute(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        startAtt = self.runner.pfa.getAttribute('peptide_start')
        startAtt.name = "peptide_start_fake"
        afw = pcssIO.AnnotationFileWriter(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            afw.writeAllOutput(self.proteins)
        self.handleTestException(e)
    
    def test_bad_get_protein_attribute(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            self.proteins[0].setStringAttribute("fake", "fakeValue")
        self.handleTestException(e)
    
    def test_read_attributes(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEquals(self.runner.pfa.getColumnSortedInputAttributes()[-1].name, "peptide_errors")

    def test_write_normal_output(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        modelColumns = pcssModels.PcssModelTableColumns(self.pcssConfig)
        self.modelTable = pcssModels.PcssModelTable(self.runner, modelColumns)
        pcssProtein = self.getProtein("76c3a409540532138c6b44bde9e4d248MDDRDENQ", self.proteins)

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

        afw = pcssIO.AnnotationFileWriter(self.runner)
        afw.writeAllOutput(self.proteins)

        annotationFileName = self.runner.pdh.getFullOutputFile("annotationOutput.txt")
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile(annotationFileName)
        newProteins = reader.getProteins()

        for newProtein in newProteins:

            oldProtein = self.getProtein(newProtein.modbaseSequenceId, self.proteins)
            self.assertTrue(oldProtein.isEqual(newProtein))
        
    def test_read_feature_error(self):
        reader = pcssIO.AnnotationFileReader(self.runner)

        reader.readAnnotationFile(self.getErrorInputFile("annotationOutputFeatureError.txt"))
        
        firstProtein = reader.getProteins()[0]
        firstPeptide = firstProtein.peptides.values()[0]
        self.assertEqual("peptide_error_no_source_model", firstPeptide.getAttributeOutputString("dssp_structure"))

    def test_no_peptides_parsed(self):
        self.runner.pcssConfig['fasta_file'] = self.getErrorInputFile("noPeptidesParsedFasta.txt")
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEqual(self.proteins[0].getAttributeOutputString("protein_errors"), "no_peptides_parsed")
        
    def test_annotation_no_sequences(self):

        self.runner.pcssConfig['fasta_file'] = self.getErrorInputFile("annotationMissingSequence.txt")
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile("testInput/svmApplicationAnnotationInput.txt")
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readProteinSequences(self.runner.pcssConfig['fasta_file'])
        self.handleTestException(e)
        

    def test_scan_peptides(self):
        self.proteins = self.spi.readInputFile(self.runner.pcssConfig['fasta_file'])
        self.assertEqual(self.proteins[0].modbaseSequenceId, "76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        self.assertEqual(self.proteins[0].uniprotId, "P62258")
        self.assertEqual(self.proteins[0].peptides.values()[0].startPosition, 2)
        self.assertEqual(self.proteins[0].peptides.values()[0].endPosition, 9)
        self.assertEqual(self.proteins[0].peptides.values()[0].sequence, "DREDLVYQ")
        self.assertEqual(len(self.proteins[0].peptides.values()), 19)
        

    def test_defined_peptides(self):
        self.runner.pcssConfig['fasta_file'] = "testInput/inputSequenceDefined.txt"
        dpi = pcssIO.DefinedPeptideImporter(self.runner)
    
        self.proteins = dpi.readInputFile(self.runner.pcssConfig['fasta_file'])
        
        self.assertEqual(self.proteins[0].modbaseSequenceId, "76c3a409540532138c6b44bde9e4d248MDDRDENQ")
        self.assertEqual(self.proteins[0].uniprotId, "P62258")
        self.assertEqual(self.proteins[0].peptides.values()[0].startPosition, 2)
        self.assertEqual(self.proteins[0].peptides.values()[0].endPosition, 9)
        self.assertEqual(self.proteins[0].peptides.values()[0].sequence, "DREDLVYQ")

    def processBadDefinedFastaTest(self, fileName):
        dpi = pcssIO.DefinedPeptideImporter(self.runner)
        
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            dpi.readInputFile(self.getErrorInputFile(fileName))
        self.handleTestException(e)

    def test_bad_defined_fasta(self):
        
        self.processBadDefinedFastaTest("badStatusFasta.txt")
        self.processBadDefinedFastaTest("noPeptidesFasta.txt")
        self.processBadDefinedFastaTest("onlyMismatchFasta.txt")
        self.processBadDefinedFastaTest("peptideMismatchFasta.txt")

    def test_no_annotation_proteins(self):
        reader = pcssIO.AnnotationFileReader(self.runner)
        with self.assertRaises(pcssErrors.PcssGlobalException) as e:
            reader.readAnnotationFile(self.getErrorInputFile("noAnnotationProteins.txt"))
        self.handleTestException(e)
        print "no annotation proteins: msg %s" % e.exception.msg

    def read_write_annotation_file(self, inputFile):
        reader = pcssIO.AnnotationFileReader(self.runner)
        reader.readAnnotationFile(inputFile)
        proteins = reader.getProteins()
        afw = pcssIO.AnnotationFileWriter(self.runner)
        afw.writeAllOutput(proteins)
        self.compareFiles(inputFile, 
                          self.runner.pdh.getFullOutputFile(self.runner.internalConfig["annotation_output_file"]), True)
        
    def test_read_write_normal_annotation_file(self):    
        self.read_write_annotation_file("testInput/svmApplicationAnnotationInput.txt")

    def test_read_write_feature_error_file(self):
        
        self.read_write_annotation_file(self.getErrorInputFile("annotationOutputFeatureError.txt"))
                    

if __name__ == '__main__':
    unittest.main()
