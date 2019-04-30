import copy
import globals

######################################################
# computeOperonArrangements
# Parameters:
# Description: Takes a list of events and computes Forward Conserved, Transposed Forward, Inverted, Inverted Transposed, and Lost regions
######################################################
def computeOperonArrangements(events):
    #Various regions of the genome
    conservedForwardRegions = []
    transposedForwardRegions = []
    invertedRegions = []
    invertedTransposedRegions = []
    lostRegions = []

    #Sort events by x coordinates
    eventsCopy = copy.deepcopy(events)
    eventsCopy.sort(key=lambda x:x.fragmentDetails1.fragmentIndex, reverse=False)

    numLossesDetected = 0 #This will be used to forgive the threshold distance between two operons in the event we encounter a lot of losses

    while len(eventsCopy) > 0:
        currEvent = eventsCopy.pop(0)

        if currEvent.score == - 1:          #This is an operon that has been lost
            lostRegions.append(currEvent)
            numLossesDetected += 1
        else:                               #This an operon on the dot plot
            consecutiveRegion = [] #Stores a consecutive region
            consecutiveRegion.append(currEvent)
            foundNeighbor = True #Keeps track of whether we found a neighboring operon

            yIncreaseCounter = 0 #Counts the number of points that are increasing in the current region
            yDecreaseCounter = 0 #Counts the number of points that are decreasing in the current region

            minMainDiagonalDistance = abs(currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex) #Calculates the distance from the main diagonal
            
            #These track whether we cross the main diagonal
            aboveMainDiagonal = False
            belowMainDiagonal = False
            if currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex < 0:
                aboveMainDiagonal = True
            if currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex > 0:
                belowMainDiagonal = True

            while foundNeighbor and len(eventsCopy) > 0: #Continue looping as long as we keep on finding a neighbor
                foundNeighbor = False #Reset the tracker

                prevEvent = consecutiveRegion[len(consecutiveRegion)-1] #Gets the last operon in the the consecutive operons list
                currEvent = eventsCopy[0] #Get the next available operon
                yDistance = abs(currEvent.fragmentDetails2.fragmentIndex - prevEvent.fragmentDetails2.fragmentIndex) #The distance on the y-axis

                if yDistance - numLossesDetected < globals.yDistanceThreshold and currEvent.score != -1: #If the y-Distance is less than the threshold and not a lost operon then add it to the consecutive region
                    consecutiveRegion.append(eventsCopy.pop(0))
                    foundNeighbor = True #Indicates we found a consecutive region
                    
                    #Tracks the minimuim distance from the main diagonal
                    currMainDiagonalDistance = abs(currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex)
                    if currMainDiagonalDistance < minMainDiagonalDistance:
                        minMainDiagonalDistance = currMainDiagonalDistance

                    #Tracks whether this point is above or below the main diagonal
                    if currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex < 0:
                        aboveMainDiagonal = True
                    if currEvent.fragmentDetails1.fragmentIndex - currEvent.fragmentDetails2.fragmentIndex > 0:
                        belowMainDiagonal = True

                    #Determines if the points are moving up or down
                    if prevEvent.fragmentDetails2.fragmentIndex < currEvent.fragmentDetails2.fragmentIndex:
                        yIncreaseCounter += 1
                    else:
                        yDecreaseCounter += 1
                elif currEvent.score == -1: #Ignore the fragments with a -1 score
                    lostRegions.append(eventsCopy.pop(0))
                    foundNeighbor = True
                    numLossesDetected +=1

            if yIncreaseCounter > yDecreaseCounter or len(consecutiveRegion) == 1:  #If the points were going up or there was only 1 point
                if minMainDiagonalDistance < globals.yDistanceThreshold:            #If the distance from the main diagonal is below the threshold
                    conservedForwardRegions.append(consecutiveRegion)
                else:                                                               #Else this is a transposition
                    transposedForwardRegions.append(consecutiveRegion)
            else:
                if aboveMainDiagonal == True and belowMainDiagonal == True:         #If the points crossed the main diagonal it's an inversion
                    invertedRegions.append(consecutiveRegion)
                else:                                                               #Else it's a transposition
                    invertedTransposedRegions.append(consecutiveRegion)

    print('Statistics for regions in the genome:')
    print('Total number of tracking events: %s' % (len(events)))
    print('Total number of forward conserved regions: %s' % len(conservedForwardRegions))
    print('Total number of forward transposed regions: %s' % len(transposedForwardRegions))
    print('Total number of inverted regions: %s' % len(invertedRegions))
    print('Total number of inverted transposed regions: %s' % len(invertedTransposedRegions))
    print('Total number of lost operons: %s' % len(lostRegions))

    #Prints indexs of the various regions computed
    for region in conservedForwardRegions:
        print('Forward Conserved Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.fragmentIndex, region[x].fragmentDetails2.fragmentIndex))
    for region in transposedForwardRegions:
        print('Forward Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.fragmentIndex, region[x].fragmentDetails2.fragmentIndex))
    for region in invertedRegions:
        print('Inverted Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.fragmentIndex, region[x].fragmentDetails2.fragmentIndex))
    for region in invertedTransposedRegions:
        print('Inverted Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.fragmentIndex, region[x].fragmentDetails2.fragmentIndex))
    for operon in lostRegions:
        print('Lost operon')
        print('%s, %s' %(operon.fragmentDetails1.fragmentIndex, operon.fragmentDetails2.fragmentIndex))

    return conservedForwardRegions, transposedForwardRegions, invertedRegions, invertedTransposedRegions, lostRegions

######################################################
# computeRegionDetails
# Parameters:
# Description: Generates a string containing all of the details about a particular region
######################################################
def computeRegionDetails(regions, description):
    temp1 = description
    temp2 = description
    
    for region in regions:
        for fragment in region:
            startPos1 = fragment.fragmentDetails1.startPositionInGenome
            startPos2 = fragment.fragmentDetails2.startPositionInGenome
            
            seq1 = copy.deepcopy(fragment.fragmentDetails1.sequence)
            seq2 = copy.deepcopy(fragment.fragmentDetails2.sequence)
            
            if fragment.fragmentDetails1.isNegativeOrientation == True:
                seq1.reverse()
            if fragment.fragmentDetails2.isNegativeOrientation == True:
                seq2.reverse()
                
            for y in range(0, len(seq1)):
                currPos = startPos1 + y
                currGene = seq1[y]
                temp1 += currGene + ' ' + str(currPos) + ', '
                
            for y in range(0, len(seq2)):
                currPos = startPos2 + y
                currGene = seq2[y]
                temp2 += currGene + ' ' + str(currPos) + ', '

    return temp1, temp2

######################################################
# determineAncestralFragmentArrangementUsingNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the sibling and the neighbor
######################################################
def determineAncestralFragmentArrangementUsingNeighbor(FCR, TR, IR, ITR, LR, NFCR, NTR, NIR, NITR, NLR):
    ancestralFragments = None
    arrangedFragments = {}
    
    #De
    for regions in FCR:
        for region in regions:
            for fragment in region:
                if fragment.fragmentDetails1.fragmentIndex == fragment.fragmentDetails2.fragmentIndex:
                    arrangedFragments[fragment.fragmentDetails1.fragmentIndex] = fragment #If the position has been preserved in both regions then keep them in same position
                
    
    
    return ancestralFragments

######################################################
# determineAncestralFragmentArrangementWithoutNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the siblings
######################################################
def determineAncestralFragmentArrangementWithoutNeighbor(FCR, TR, IR, ITR, LR):
    ancestralFragments = None
    
    return ancestralFragments