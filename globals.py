##################################################
###Global variables used in the script############
##################################################
def initialize():
    global trackingId
    trackingId = 0

    global codonCost
    codonCost = 0.5

    global deletionCost
    deletionCost = 1

    global substitutionCost
    substitutionCost = 1
    
    global sizeDuplications
    sizeDuplications = {}
    
    global sizeDeletions
    sizeDeletions = {}
    
    global localSizeDuplications
    localSizeDuplications = {}
    
    global localSizeDeletions
    localSizeDeletions = {}
    
    global yDistanceThreshold
    yDistanceThreshold = 3