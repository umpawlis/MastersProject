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
    events = []
    deletionSizes = {}
    duplicationSizes = {}

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
    def getEvents(self):
        return self.events

    def setEvents(self, events):
        self.events = events
    def setDeletionSizes(self, deletionSizes):
        self.deletionSizes = deletionSizes
    def setDuplicationSizes(self, duplicationSizes):
        self.duplicationSizes = duplicationSizes