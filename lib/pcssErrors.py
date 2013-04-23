import sys
import os
from Bio import SeqIO
import itertools
import pcssTools
import logging
log = logging.getLogger("pcssPeptide")


class PcssShutilError(Exception):
    def __init__(self, ioError, function, args):
        self.ioError = ioError
        self.function = function.__name__
        self.args = args

class PcssGlobalException(Exception):    
    def __init__(self,  msg):
        self.msg = msg

class ProteinException(Exception):
    def __init__(self, msg):
        self.msg = msg

    def getErrorCodePrefix(self):
        return pcssTools.getProteinErrorCodePrefix()

    def setProtein(self, pcssProtein):
        self.pcssProtein = pcssProtein

class ErrorInfo:
    def __init__(self, errorFile):
        reader = pcssTools.PcssFileReader(errorFile)
        lines = reader.getLines()
        self.errorType = lines[0]
        errorLine = lines[1]
        self.msg = errorLine


class ErrorExistsException(Exception):
    def __init__(self, msg, fileName):
        self.msg = msg
        self.fileName = fileName

class InternalException(Exception):
    def __init__(self, msg):
        self.msg = msg

class PeptideException(ProteinException):
    def __init__(self,  msg):
        self.msg = msg


    def getErrorCodePrefix(self):
        return pcssTools.getPeptideErrorCodePrefix()

    
    def setPeptide(self, pcssPeptide):
        self.pcssPeptide = pcssPeptide


    

class PsipredException(PeptideException):
    def __init__(self, msg):
        self.msg = msg

    def getAffectedAttributeNames(self):
        return ["psipred_string_feature", "psipred_score_feature"]

class DisopredException(PeptideException):
    def __init__(self, msg):
        self.msg = msg

    def getAffectedAttributeNames(self):
        return ["disopred_string_feature", "disopred_score_feature"]



class DisopredPeptideNotFoundException(DisopredException):
    def __init__(self, msg):
        self.msg = msg
        self.code = self.getErrorCodePrefix() + "disopred_peptide_not_found"

class PsipredPeptideNotFoundException(PsipredException):
    def __init__(self, msg):
        self.msg = msg
        self.code = self.getErrorCodePrefix() + "psipred_peptide_not_found"

class DisopredMismatchException(DisopredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "disopred_protein_mismatch"
        self.msg = msg

class PsipredMismatchException(PsipredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "psipred_protein_mismatch"
        self.msg = msg

class DisopredBadCallException(DisopredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "disopred_bad_call"
        self.msg = msg

class DisopredCommandException(DisopredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "disopred_bad_command"
        self.msg = msg

class DisopredBadLineException(DisopredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "disopred_bad_line"
        self.msg = msg

class PsipredBadLineException(PsipredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "psipred_bad_line"
        self.msg = msg

class PsipredBadCallException(PsipredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "psipred_bad_call"
        self.msg = msg

class PsipredCommandException(PsipredException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "psipred_bad_command"
        self.msg = msg

class StructureException(PeptideException):
    def __init__(self, msg):
        self.msg = msg

    def setModel(self, model):
        self.model = model

    def getAffectedAttributeNames(self):
        return ["dssp_structure", "dssp_accessibility"]

class NoSourceModelException(StructureException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "no_source_model"
        self.msg = msg
        

class DsspException(StructureException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "dssp_error"
        self.msg = msg


class DsspMismatchException(StructureException):
    def __init__(self, msg):
        self.code = self.getErrorCodePrefix() + "dssp_mismatch"
        self.msg = msg
    
    
