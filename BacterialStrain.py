######################################################
# BacterialStrain
# Parameters:
# Description: Stores information about the specific bacterial strain
######################################################
class BacterialStrain(object):

    #Class constructor
    def __init__(self, name, genomeFragments):
        self.name = name
        self.genomeFragments = genomeFragments