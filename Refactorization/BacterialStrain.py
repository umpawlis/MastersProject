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
        
    def addCodonMismatchDetails(self, codonMismatchDetails):
        self.codonMismatchDetails += codonMismatchDetails
        
    def addSubstitutionDetails(self, substitutionDetails):
        self.substitutionDetails += substitutionDetails