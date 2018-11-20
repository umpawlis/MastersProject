######################################################
# Tracking Event
# Parameters:
# Description: Stores information about a pair of orthologous operons
######################################################
class TrackingEvent(object):
    trackingEventId = 0
    score = -1
    genome1Name = ""
    genome2Name = ""
    genome1Operon = ""
    genome2Operon = ""
    genome1OperonIndex = -1
    genome2OperonIndex = -1
    ancestralOperon = ""
    technique = ""
    numLosses = 0
    operonEvents = None

    def __init__(self, trackingEventId, score, genome1Name, genome2Name, genome1Operon, genome2Operon, genome1OperonIndex, genome2OperonIndex, ancestralOperon, technique):
        self.trackingEventId = trackingEventId
        self.score = score
        self.genome1Name = genome1Name
        self.genome2Name = genome2Name
        self.genome1Operon = genome1Operon
        self.genome2Operon = genome2Operon
        self.genome1OperonIndex = genome1OperonIndex
        self.genome2OperonIndex = genome2OperonIndex
        self.ancestralOperon = ancestralOperon
        self.technique = technique
        self.lostEventIds = []

    def printTrackingEvent(self):
        if self.operonEvents != None:
            operonEventsString = self.operonEvents.outputOperonEvent()
            print("{ Tracking Event Id: %s, \nAlignment Score: %s, \nGenome 1: %s, \nGenome 2: %s, \nGenome 1 Operon: %s, \nGenome 2 Operon: %s, \nGenome 1 Operon Index: %s, \nGenome 2 Operon Index: %s, \nAncestral Operon: %s, \nTechnique: %s, \nLost Event Ids: %s, \nOperon Events: %s}"%(self.trackingEventId, self.score, self.genome1Name, self.genome2Name, self.genome1Operon, self.genome2Operon, self.genome1OperonIndex, self.genome2OperonIndex, self.ancestralOperon, self.technique, self.lostEventIds, operonEventsString))
        else:
            print("{ Tracking Event Id: %s, \nAlignment Score: %s, \nGenome 1: %s, \nGenome 2: %s, \nGenome 1 Operon: %s, \nGenome 2 Operon: %s, \nGenome 1 Operon Index: %s, \nGenome 2 Operon Index: %s, \nAncestral Operon: %s, \nTechnique: %s, \nLost Event Ids: %s}"%(self.trackingEventId, self.score, self.genome1Name, self.genome2Name, self.genome1Operon, self.genome2Operon, self.genome1OperonIndex, self.genome2OperonIndex, self.ancestralOperon, self.technique, self.lostEventIds))
    #####################################
    #Getters
    #####################################
    def getTrackingEventId(self):
        return self.trackingEventId
    def getScore(self):
        return self.score
    def getGenome1Name(self):
        return self.genome1Name
    def getGenome2Name(self):
        return self.genome2Name
    def getGenome1Operon(self):
        return self.genome1Operon
    def getGenome2Operon(self):
        return self.genome2Operon
    def getGenome1OperonIndex(self):
        return self.genome1OperonIndex
    def getGenome2OperonIndex(self):
        return self.genome2OperonIndex
    def getAncestralOperon(self):
        return self.ancestralOperon
    def getTechnique(self):
        return self.technique
    def getLostEventIds(self):
        return self.lostEventIds
    def getOperonEvents(self):
        return self.operonEvents

    #####################################
    #Setters
    #####################################
    def setTrackingEventId(self, trackingEventId):
        self.trackingEventId = trackingEventId
    def setScore(self, score):
        self.score = score
    def setGenome1Name(self, genome1Name):
        self.genome1Name = genome1Name
    def setGenome2Name(self, genome2Name):
        self.genome2Name = genome2Name
    def setGenome1Operon(self, genome1Operon):
        self.genome1Operon = genome1Operon
    def setGenome2Operon(self, genome2Operon):
        self.genome2Operon = genome2Operon
    def setGenome1OperonIndex(self, genome1OperonIndex):
        self.genome1OperonIndex = genome1OperonIndex
    def setGenome2OperonIndex(self, genome2OperonIndex):
        self.genome2OperonIndex = genome2OperonIndex
    def setAncestralOperon(self, ancestralOperon):
        self.ancestralOperon = ancestralOperon
    def setTechnique(self, technique):
        self.technique = technique
    def setLostEventIds(self, lostEventIds):
        self.lostEventIds = lostEventIds
    def setOperonEvents(self, operonEvents):
        self.operonEvents = operonEvents