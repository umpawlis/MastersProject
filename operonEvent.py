######################################################
# Operon Events
# Parameters:
# Description: Stores information about the events that have occured in the operon
######################################################
class OperonEvents(object):
    numMatches = 0
    numCodonMismatches = 0
    numMismatches = 0
    numSubstitutions = 0
    operon1 = None
    operon2 = None
    matrix = None
    operon1GeneLosses = 0
    operon2GeneLosses = 0
    operon1GeneDuplicates = 0
    operon2GeneDuplicates = 0
    operon1Gaps = []
    operon2Gaps = []
    duplicateSizesOp1 = []
    duplicateSizesOp2 = []
    alignedGenesInOperon1 = []
    alignedGenesInOperon2 = []
    lossesDueToSlidingWindowMethodOperon1 = 0
    lossesDueToSlidingWindowMethodOperon2 = 0
    duplicationsDueToSlidingWindowMethodOperon1 = 0
    duplicationsDueToSlidingWindowMethodOperon2 = 0
    operon1GapIndexes = []
    operon2GapIndexes = []
    substitutionDict = {}
    operon1IndexToAlignment2Index = {}
    operon1Reversed = False
    operon2Reversed = False

    def __init__(self, numMatches, numCodonMismatches, numMismatches, numSubstitutions, operon1, operon2, matrix):
        self.numMatches = numMatches
        self.numCodonMismatches = numCodonMismatches
        self.numMismatches = numMismatches
        self.numSubstitutions = numSubstitutions
        self.operon1 = operon1
        self.operon2 = operon2
        self.matrix = matrix

    def toStringOperonEvents(self):
        return "(Num Matches = %s,\nNum Codon Mismatches = %s,\nNum Mismatches = %s,\nNum Substitutions = %s,\nOperon 1 = %s,\nOperon 2 = %s,\nOperon 1 Gene Losses: %s,\nOperon 2 Gene Losses: %s,\nOperon 1 Gene Duplications: %s,\nOperon 2 Gene Duplications: %s)" % (self.numMatches, self.numCodonMismatches, self.numMismatches, self.numSubstitutions, self.operon1, self.operon2, self.operon1GeneLosses, self.operon2GeneLosses, self.operon1GeneDuplicates, self.operon2GeneDuplicates)

    def toStringAlignmentResults(self):
        return "Results of Alignment: (Num Matches = %s,\nNum Codon Mismatches = %s,\nNum Mismatches = %s,\nNum Substitutions = %s,\nOperon 1 = %s,\nOperon 2 = %s)" % (self.numMatches, self.numCodonMismatches, self.numMismatches, self.numSubstitutions, self.operon1, self.operon2)

    def outputOperonEvent(self):
        return ("(Operons Compared = %s, %s\nNumber of Matches = %s\nNumber of Codon Mismatches = %s\nNumber of Mismatches = %s\nNumber of Substitutions = %s\nExtra Genes In Operon 1 = %s\nExtra Genes in Operon 2 = %s\nAlignment in Operon 1 = %s\nAlignment in Operon 2 = %s)" % (self.operon1, self.operon2, self.numMatches, self.numCodonMismatches, self.numMismatches, self.numSubstitutions, self.operon1Gaps, self.operon2Gaps, self.alignedGenesInOperon1, self.alignedGenesInOperon2) )
        
    #####Getters#####
    def getNumMatches(self):
        return self.numMatches
    def getNumCodonMismatches(self):
        return self.numCodonMismatches
    def getNumMismatches(self):
        return self.numMismatches
    def getNumSubstitutions(self):
        return self.numSubstitutions
    def getOperon1(self):
        return self.operon1
    def getOperon2(self):
        return self.operon2
    def getMatrix(self):
        return self.matrix
    def getOperon1GeneLosses(self):
        return self.operon1GeneLosses
    def getOperon2GeneLosses(self):
        return self.operon2GeneLosses
    def getOperon1GeneDuplicates(self):
        return self.operon1GeneDuplicates
    def getOperon2GeneDuplicates(self):
        return self.operon2GeneDuplicates
    def getOperon1Gaps(self):
        return self.operon1Gaps
    def getOperon2Gaps(self):
        return self.operon2Gaps
    def getDuplicateSizesOp1(self):
        return self.duplicateSizesOp1
    def getDuplicateSizesOp2(self):
        return self.duplicateSizesOp2
    def getAlignedGenesInOperon1(self):
        return self.alignedGenesInOperon1
    def getAlignedGenesInOperon2(self):
        return self.alignedGenesInOperon2
    def getLossesDueToSlidingWindowMethodOperon1(self):
        return self.lossesDueToSlidingWindowMethodOperon1
    def getLossesDueToSlidingWindowMethodOperon2(self):
        return self.lossesDueToSlidingWindowMethodOperon2
    def getDuplicationsDueToSlidingWindowMethodOperon1(self):
        return self.duplicationsDueToSlidingWindowMethodOperon1
    def getDuplicationsDueToSlidingWindowMethodOperon2(self):
        return self.duplicationsDueToSlidingWindowMethodOperon2
    def getOperon1GapIndexes(self):
        return self.operon1GapIndexes
    def getOperon2GapIndexes(self):
        return self.operon2GapIndexes
    def getSubstitutionDict(self):
        return self.substitutionDict
    def getOperon1IndexToAlignment2Index(self):
        return self.operon1IndexToAlignment2Index
    def isOperon1Reversed(self):
        return self.operon1Reversed
    def isOperon2Reversed(self):
        return self.operon2Reversed
    
    #####Setters#####
    def setNumMatches(self, numMatches):
        self.numMatches = numMatches
    def setNumCodonMismatches(self, numCodonMismatches):
        self.numCodonMismatches = numCodonMismatches
    def setNumMismatches(self, numMismatches):
        self.numMismatches = numMismatches
    def setNumSubstitutions(self, numSubstitutions):
        self.numSubstitutions = numSubstitutions
    def setOperon1(self, operon1):
        self.operon1 = operon1
    def setOperon2(self, operon2):
        self.operon2 = operon2
    def setMatrix(self, matrix):
        self.matrix = matrix
    def setOperon1GeneLosses(self, operon1GeneLosses):
        self.operon1GeneLosses = operon1GeneLosses
    def setOperon2GeneLosses(self, operon2GeneLosses):
        self.operon2GeneLosses = operon2GeneLosses
    def setOperon1GeneDuplicates(self, operon1GeneDuplicates):
        self.operon1GeneDuplicates = operon1GeneDuplicates
    def setOperon2GeneDuplicates(self, operon2GeneDuplicates):
        self.operon2GeneDuplicates = operon2GeneDuplicates
    def setOperon1Gaps(self, operon1Gaps):
        self.operon1Gaps = operon1Gaps
    def setOperon2Gaps(self, operon2Gaps):
        self.operon2Gaps = operon2Gaps
    def setDuplicateSizesOp1(self, duplicateSizesOp1):
        self.duplicateSizesOp1 = duplicateSizesOp1
    def setDuplicateSizesOp2(self, duplicateSizesOp2):
        self.duplicateSizesOp2 = duplicateSizesOp2
    def setAlignedGenesInOperon1(self, alignedGenesInOperon1):
        self.alignedGenesInOperon1 = alignedGenesInOperon1
    def setAlignedGenesInOperon2(self, alignedGenesInOperon2):
        self.alignedGenesInOperon2 = alignedGenesInOperon2
    def setLossesDueToSlidingWindowMethodOperon1(self, lossesDueToSlidingWindowMethodOperon1):
        self.lossesDueToSlidingWindowMethodOperon1 = lossesDueToSlidingWindowMethodOperon1
    def setLossesDueToSlidingWindowMethodOperon2(self, lossesDueToSlidingWindowMethodOperon2):
        self.lossesDueToSlidingWindowMethodOperon2 = lossesDueToSlidingWindowMethodOperon2
    def setDuplicationsDueToSlidingWindowMethodOperon1(self, duplicationsDueToSlidingWindowMethodOperon1):
        self.duplicationsDueToSlidingWindowMethodOperon1 = duplicationsDueToSlidingWindowMethodOperon1
    def setDuplicationsDueToSlidingWindowMethodOperon2(self, duplicationsDueToSlidingWindowMethodOperon2):
        self.duplicationsDueToSlidingWindowMethodOperon2 = duplicationsDueToSlidingWindowMethodOperon2
    def setOperon1GapIndexes(self, operon1GapIndexes):
        self.operon1GapIndexes = operon1GapIndexes
    def setOperon2GapIndexes(self, operon2GapIndexes):
        self.operon2GapIndexes = operon2GapIndexes
    def setSubstitutionDict(self, substitutionDict):
        self.substitutionDict = substitutionDict
    def setOperon1IndexToAlignment2Index(self, operon1IndexToAlignment2Index):
        self.operon1IndexToAlignment2Index = operon1IndexToAlignment2Index
    def setOperon1Reversed(self, operon1Reversed):
        self.operon1Reversed = operon1Reversed
    def setOperon2Reversed(self, operon2Reversed):
        self.operon2Reversed = operon2Reversed