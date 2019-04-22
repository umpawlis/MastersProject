######################################################
# BacterialStrain
# Parameters:
# Description: Stores information about the specific bacterial strain
######################################################
class BacterialStrain(object):

    #Class constructor
    def __init__(self, name, genomeFragments):
        self.name = name
        self.genomeFragments = genomeFragments
        self.codonMismatchDetails = 'Codon Mismatch:'
        self.substitutionDetails = 'Substitution:'
        self.duplicationCounts = {}
        self.duplicationDetails = 'Duplication:'
        self.deletionCounts = {}
        self.deletionDetails = 'Deletion:'
        
    def addCodonMismatchDetails(self, codonMismatchDetails):
        self.codonMismatchDetails += codonMismatchDetails
        
    def addSubstitutionDetails(self, substitutionDetails):
        self.substitutionDetails += substitutionDetails