class Event(object):

    def __init__(self, trackingEventId):
        self.trackingEventId = trackingEventId

    def setScore(self, score):
        self.score = score

    def isOriginallyNegativeOrientationOp1(self, originallyNegativeOrientationOp1):
        self.originallyNegativeOrientationOp1 = originallyNegativeOrientationOp1

    def isOriginallyNegativeOrientationOp2(self, originallyNegativeOrientationOp2):
        self.originallyNegativeOrientationOp2 = originallyNegativeOrientationOp2

    def isReversedOp1(self, reversedOp1):
        self.reversedOp1 = reversedOp1

    def isReversedOp2(self, reversedOp2):
        self.reversedOp2 = reversedOp2

    def setGenome1Name(self, genome1Name):
        self.genome1Name = genome1Name

    def setGenome2Name(self, genome2Name):
        self.genome2Name = genome2Name

    def setGenome1Operon(self, genome1Operon):
        self.genome1Operon = genome1Operon

    def setGenome2Operon(self, genome2Operon):
        self.genome2Operon = genome2Operon

    def setNumMatches(self, numMatches):
        self.numMatches = numMatches

    def setNumCodonMismatches(self, numCodonMismatches):
        self.numCodonMismatches = numCodonMismatches

    def setNumGeneMismatches(self, numGeneMismatches):
        self.numGeneMismatches = numGeneMismatches

    def setNumSubstitutions(self, numSubstitutions):
        self.numSubstitutions = numSubstitutions

    def setOperon1Alignment(self, operon1Alignment):
        self.operon1Alignment = operon1Alignment

    def setOperon2Alignment(self, operon2Alignment):
        self.operon2Alignment = operon2Alignment

    def setOperon1Gaps(self, operon1Gaps):
        self.operon1Gaps = operon1Gaps

    def setOperon2Gaps(self, operon2Gaps):
        self.operon2Gaps = operon2Gaps

    def setOperon1GapIndexes(self, operon1GapIndexes):
        self.operon1GapIndexes = operon1GapIndexes

    def setOperon2GapIndexes(self, operon2GapIndexes):
        self.operon2GapIndexes = operon2GapIndexes

    def setOperon1Index(self, operon1Index):
        self.operon1Index = operon1Index

    def setOperon2Index(self, operon2Index):
        self.operon2Index = operon2Index

    def setTechnique(self, technique):
        self.technique = technique
        
    def setAncestralOperonGeneSequence(self, ancestralOperonGeneSequence):
        self.ancestralOperonGeneSequence = ancestralOperonGeneSequence
        
    def setSizeDuplications(self, sizeDuplications):
        self.sizeDuplications = sizeDuplications
    
    def setSizeDeletions(self, sizeDeletions):
        self.sizeDeletions = sizeDeletions

    def printEvent(self):
        print('Tracking Id: %s\nTwo strains compared: %s, %s\nScore: %s\nOperon 1: %s\nOperon 2: %s\nOperon 1 position: %s\nOperon 2 position: %s\nAncestral Operon Sequence: %s' % (self.trackingEventId, self.genome1Name, self.genome2Name, self.score, self.genome1Operon, self.genome2Operon, self.operon1Index, self.operon2Index, self.ancestralOperonGeneSequence))