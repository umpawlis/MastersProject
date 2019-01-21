
######################################################
# Lineage Data
# Parameters:
# Description: Stores data about a branch
######################################################
class LineageData(object):
    name = ""
    geneDuplications = 0
    geneLosses = 0
    operonDuplications = {}
    operonLosses = {}
    lineageCost = 0
    
    def __init__(self, name, geneDuplications, geneLosses, operonDuplications, operonLosses, lineageCost):
        self.name = name
        self.geneDuplications = geneDuplications
        self.geneLosses = geneLosses
        self.operonDuplications = operonDuplications
        self.operonLosses = operonLosses
        self.lineageCost = lineageCost
        
    #Getters
    def getName(self):
        return self.name
    def getGeneDuplications(self):
        return self.geneDuplications
    def getGeneLosses(self):
        return self.geneLosses
    def getOperonDuplications(self):
        return self.operonDuplications
    def getOperonLosses(self):
        return self.operonLosses
    def getLineageCost(self):
        return self.lineageCost
    
    #Setters
    def setName(self, name):
        self.name = name
    def setGeneDuplications(self, geneDuplications):
        self.geneDuplications = geneDuplications
    def setGeneLosses(self, geneLosses):
        self.geneLosses = geneLosses
    def setOperonDuplications(self, operonDuplications):
        self.operonDuplications = operonDuplications
    def setOperonLosses(self, operonLosses):
        self.operonLosses = operonLosses
    def setLineageCost(self, lineageCost):
        self.lineageCost = lineageCost
        
    #ToString