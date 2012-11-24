import sys
import os
import re
import subprocess
import configobj
import pcssTools
import pcssPeptide


configFile = sys.argv[1]

pcssConfig = configobj.ConfigObj(configFile)

runner = pcssTools.PcssRunner(pcssConfig)

spi = pcssPeptide.ScanPeptideImporter()
spi.readPeptides(pcssConfig['fasta_file'])
