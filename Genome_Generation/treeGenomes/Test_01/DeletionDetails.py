######################################################
# DeletionDetails
# Parameters:
# Description: Stores information about deleted genes that we may want to mark as duplications later on
######################################################
class DeletionDetails(object):

    #Class constructor
    def __init__(self, ancestralGene, ancestralPosition, strainFragmentId, strain, originalGene, originalPosition, ancestralOperon, otherStrain):
        self.ancestralGene = ancestralGene
        self.ancestralPosition = ancestralPosition
        self.strainFragmentId = strainFragmentId
        self.strain = strain
        self.originalGene = originalGene
        self.originalPosition = originalPosition
        self.ancestralOperon = ancestralOperon
        self.ancestralFragmentId = -1
        self.geneRemoved = False
        self.otherStrain = otherStrain
        
    def toString(self):
        return "Ancestral Operon: %s\nAncestral Gene: %s\nAncestral Position: %s\nFragment Id: %s\nOriginal Gene: %s\nOriginal Position: %s\nOriginal Strain: %s" %(self.ancestralOperon, self.ancestralGene, self.ancestralPosition, self.strainFragmentId, self.originalGene, self.originalPosition, self.strain.name)