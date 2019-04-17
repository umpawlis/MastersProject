######################################################
# BacterialStrain
# Parameters:
# Description: Stores information about the specific bacterial strain
######################################################
class BacterialStrain(object):

    #Class constructor
    def __init__(self, name, genomeFragments, formattedSequence, sequenceConversion):
        self.name = name
        self.genomeFragments = genomeFragments
        self.formattedSequence = formattedSequence
        self.sequenceConversion = sequenceConversion
        self.deletionSizes = {}
        self.duplicationSizes = {}
        self.deletionInfo = "Deletion:"
        self.duplicationInfo = "Duplication:"
        self.codonMismatchDetails = "Codon Mismatch:"
        self.numCodonMismatches = 0