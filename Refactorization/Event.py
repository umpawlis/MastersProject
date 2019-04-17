class Event(object):

    def __init__(self, trackingEventId):
        self.trackingEventId = trackingEventId

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
        