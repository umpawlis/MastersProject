######################################################
# LineageSummary
# Parameters:
# Description: Stores the cost of a lineage
######################################################
class LineageSummary(object):
    #Class constructor
    def __init__(self, name):
        self.name = name
        self.totalCodonMismatches = 0
        self.totalSubstitutions = 0