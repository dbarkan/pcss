import sys
import os
import re
import pcssTools


if (len(sys.argv) < 4):
    print "Usage: python getTwoLetterFasta.py <sourceFastaName> <twoLetterCode> <outputFastaName>\n"
    print "Simple script to grab all fasta entries starting with modbase sequence id <twoLetterCode> from <sourceFastaName>"
    print "and output to <outputFastaName>"
    sys.exit()

sourceFastaName = sys.argv[1]
twoLetterCode = sys.argv[2]
outputFastaName = sys.argv[3]

grabber = pcssTools.TwoLetterFastaGrabber(sourceFastaName, twoLetterCode, outputFastaName)
