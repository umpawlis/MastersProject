######################################################
# GenomeFragment
# Parameters:
# Description: Stores information about either an operon or a singleton gene
######################################################
class GenomeFragment(object):

    #Class constructor
    def __init__(self, fragmentIndex, originalSequence, sequence, startPositionInGenome, description, isNegativeOrientation):
        self.fragmentIndex = fragmentIndex
        self.originalSequence = originalSequence
        self.sequence = sequence
        self.startPositionInGenome = startPositionInGenome
        self.description = description
        self.isNegativeOrientation = isNegativeOrientation
        self.deletionDetailsList = []
        self.isDuplicate = False

    def setPoint(self, point):
        self.point = point