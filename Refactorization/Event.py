class Event(object):

    def __init__(self, trackingEventId):
        self.trackingEventId = trackingEventId
        self.numCodonMismatches = 0
        self.numSubstitutions = 0

    def setScore(self, score):
        self.score = score

    def setFragmentDetails1(self, fragmentDetails1):
        self.fragmentDetails1 = fragmentDetails1

    def setFragmentDetails2(self, fragmentDetails2):
        self.fragmentDetails2 = fragmentDetails2

    def setGenome1Name(self, genome1Name):
        self.genome1Name = genome1Name

    def setGenome2Name(self, genome2Name):
        self.genome2Name = genome2Name

    def setTechnique(self, technique):
        self.technique = technique

    def setOperon1Alignment(self, operon1Alignment):
        self.operon1Alignment = operon1Alignment

    def setOperon2Alignment(self, operon2Alignment):
        self.operon2Alignment = operon2Alignment

    def setNumMatches(self, numMatches):
        self.numMatches = numMatches

    def setNumMismatches(self, numMismatches):
        self.numMismatches = numMismatches

    def setNumCodonMismatches(self, numCodonMismatches):
        self.numCodonMismatches = numCodonMismatches

    def setCodonMismatchIndexesStrain1(self, codonMismatchIndexesStrain1):
        self.codonMismatchIndexesStrain1 = codonMismatchIndexesStrain1

    def setCodonMismatchIndexesStrain2(self, codonMismatchIndexesStrain2):
        self.codonMismatchIndexesStrain2 = codonMismatchIndexesStrain2

    def setCodonMismatchGenesStrain1(self, codonMismatchGenesStrain1):
        self.codonMismatchGenesStrain1 = codonMismatchGenesStrain1

    def setCodonMismatchGenesStrain2(self, codonMismatchGenesStrain2):
        self.codonMismatchGenesStrain2 = codonMismatchGenesStrain2

    def setNumSubstitutions(self, numSubstitutions):
        self.numSubstitutions = numSubstitutions

    def setSubstitutionIndexesStrain1(self, substitutionIndexesStrain1):
        self.substitutionIndexesStrain1 = substitutionIndexesStrain1

    def setSubstitutionIndexesStrain2(self, substitutionIndexesStrain2):
        self.substitutionIndexesStrain2 = substitutionIndexesStrain2

    def setSubstitutionGenesStrain1(self, substitutionGenesStrain1):
        self.substitutionGenesStrain1 = substitutionGenesStrain1

    def setSubstitutionGenesStrain2(self, substitutionGenesStrain2):
        self.substitutionGenesStrain2 = substitutionGenesStrain2
        
    def setOperon1Gaps(self, operon1Gaps):
        self.operon1Gaps = operon1Gaps

    def setOperon2Gaps(self, operon2Gaps):
        self.operon2Gaps = operon2Gaps
        
    def setOperon1GapIndexes(self, operon1GapIndexes):
        self.operon1GapIndexes = operon1GapIndexes

    def setOperon2GapIndexes(self, operon2GapIndexes):
        self.operon2GapIndexes = operon2GapIndexes
        
    def toString(self):
        return 'Strains: %s, %s\nScore: %s\nOperon 1: %s\nOperon 2:%s\n' % (self.genome1Name, self.genome2Name, self.score, self.fragmentDetails1.originalSequence, self.fragmentDetails2.originalSequence)