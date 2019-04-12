######################################################
# GenomeFragment
# Parameters:
# Description: Stores information about either an operon or a singleton gene
######################################################
class GenomeFragment(object):

    #Class constructor
    def __init__(self, originalSequence, sequence, startPositionInGenome, description, isNegativeOrientation):
        self.originalSequence = originalSequence
        self.sequence = sequence
        self.startPositionInGenome = startPositionInGenome
        self.description = description
        self.isNegativeOrientation = isNegativeOrientation