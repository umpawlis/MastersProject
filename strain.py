######################################################
# Strain
# Parameters:
# Description: Stores information about a Strain from the newick tree
######################################################
class Strain(object):
    name = ""
    sequence = []
    genes = []
    descendants = []
    operonPositions = []
    singletonDict = {}
    trackingEvents = []

    #Class constructor
    def __init__(self, name, sequence, descendants, operonPositions, singletonDict):
        self.name = name
        self.sequence = sequence
        self.descendants = descendants
        self.operonPositions = operonPositions
        self.singletonDict = singletonDict

    #Prints the Strain content
    def printStrain(self):
        print("{ Name: %s, Sequence: %s, Descendants: %s }" %(self.name, self.sequence, self.descendants))

    def getName(self):
        return self.name
    def getSequence(self):
        return self.sequence
    def getDescendants(self):
        return self.descendants
    def getGenes(self):
        return self.genes
    def setGenes(self, genes):
        self.genes = genes
    def getOperonPositions(self):
        return self.operonPositions
    def getSingletonDict(self):
        return self.singletonDict
    def getTrackingEvents(self):
        return self.trackingEvents

    def setTrackingEvents(self, trackingEvents):
        self.trackingEvents = trackingEvents