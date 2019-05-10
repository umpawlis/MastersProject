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

    global yDistanceThreshold
    yDistanceThreshold = 3

    global ancestralCounter
    ancestralCounter = 0
    
    global inversionCounter
    inversionCounter = 0

    global transposedCounter
    transposedCounter = 0

    global invertedTransposedCounter
    invertedTransposedCounter = 0

    global deletionSizeCounter
    deletionSizeCounter = {}

    global duplicationSizeCounter
    duplicationSizeCounter = {}