import copy
import globals
from GenomeFragment import GenomeFragment

######################################################
# determineRegions
# Parameters:
# Description: Takes a list of genome fragments and computes Forward Conserved, Transposed Forward, Inverted, Inverted Transposed, and Lost regions
######################################################
def determineRegions(fragments):
    #Various regions of the genome
    conservedForwardRegions = []
    transposedForwardRegions = []
    invertedRegions = []
    invertedTransposedRegions = []

    #Sort the fragments
    fragmentsCopy = copy.deepcopy(fragments)
    fragmentsCopy.sort(key=lambda x:x.fragmentDetails1.point, reverse=False)

    while len(fragmentsCopy) > 0:
        currFragment = fragmentsCopy.pop(0)

        #Remember there's no lost fragments at this point
        consecutiveRegion = [] #Stores a consecutive region
        consecutiveRegion.append(currFragment)
        foundNeighbor = True #Keeps track of whether we found a neighboring operon

        yIncreaseCounter = 0 #Counts the number of points that are increasing in the current region
        yDecreaseCounter = 0 #Counts the number of points that are decreasing in the current region

        minMainDiagonalDistance = abs(currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point) #Calculates the distance from the main diagonal

        #These track whether we cross the main diagonal
        aboveMainDiagonal = False
        belowMainDiagonal = False
        if currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point < 0:
            aboveMainDiagonal = True
        if currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point > 0:
            belowMainDiagonal = True

        while foundNeighbor and len(fragmentsCopy) > 0: #Continue looping as long as we keep on finding a neighbor
            foundNeighbor = False #Reset the tracker

            prevFragment = consecutiveRegion[len(consecutiveRegion)-1] #Gets the last operon in the the consecutive operons list
            currFragment = fragmentsCopy[0] #Get the next available operon
            yDistance = abs(currFragment.fragmentDetails2.point - prevFragment.fragmentDetails2.point) #The distance on the y-axis

            if yDistance < globals.yDistanceThreshold: #If the y-Distance is less than the threshold add it to the consecutive region
                consecutiveRegion.append(fragmentsCopy.pop(0))
                foundNeighbor = True #Indicates we found a consecutive region

                #Tracks the minimuim distance from the main diagonal
                currMainDiagonalDistance = abs(currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point)
                if currMainDiagonalDistance < minMainDiagonalDistance:
                    minMainDiagonalDistance = currMainDiagonalDistance

                #Tracks whether this point is above or below the main diagonal
                if currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point < 0:
                    aboveMainDiagonal = True
                if currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point > 0:
                    belowMainDiagonal = True

                #Determines if the points are moving up or down
                if prevFragment.fragmentDetails2.point < currFragment.fragmentDetails2.point:
                    yIncreaseCounter += 1
                else:
                    yDecreaseCounter += 1
        #Add the region to the appropriate array
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
    print('Total number of tracking points on the graph: %s' % (len(fragments)))
    print('Total number of forward conserved regions: %s' % len(conservedForwardRegions))
    print('Total number of forward transposed regions: %s' % len(transposedForwardRegions))
    print('Total number of inverted regions: %s' % len(invertedRegions))
    print('Total number of inverted transposed regions: %s' % len(invertedTransposedRegions))

    #Prints indexs of the various regions computed
    for region in conservedForwardRegions:
        print('Forward Conserved Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point))
    for region in transposedForwardRegions:
        print('Forward Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point))
    for region in invertedRegions:
        print('Inverted Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point))
    for region in invertedTransposedRegions:
        print('Inverted Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point))

    return conservedForwardRegions, transposedForwardRegions, invertedRegions, invertedTransposedRegions

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
            temp1 += ';'

            for y in range(0, len(seq2)):
                currPos = startPos2 + y
                currGene = seq2[y]
                temp2 += currGene + ' ' + str(currPos) + ', '
            temp2 += ';'

    return temp1, temp2

######################################################
# determineAncestralFragmentArrangementUsingNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the siblings and the neighbor
######################################################
def determineAncestralFragmentArrangementUsingNeighbor(FCR, TR, IR, ITR, LR, NFCR, NTR, NIR, NITR, NLR):
    #Initialize
    ancestralFragments = []
    arrangedFragments = {}

    #Insert the lost operons at the positions they were lost in
    for fragment in LR:
        index = fragment.fragmentDetails1.fragmentIndex

        if index in arrangedFragments: #We already have a fragment at this position so just append it
            arrangedFragments[index].append(fragment)
        else: #We don't have a fragment so initialize the position
            arrangedFragments[index] = []
            arrangedFragments[index].append(fragment)

    #Insert the forward conserved regions into the new genome
    for region in FCR:
        for x in range(0, len(region)):
            fragment = region[x]
            index1 = fragment.fragmentDetails1.fragmentIndex
            index2 = fragment.fragmentDetails2.fragmentIndex

            if (index1 == index2): #The positions of the fragments has been conserved in both genomes
                if index1 in arrangedFragments:
                    arrangedFragments[index1].append(fragment)
                else:
                    arrangedFragments[index1] = []
                    arrangedFragments[index1].append(fragment)
            else: #The positions have not been conserved, means there's some operons lost between the two genomes. Insert into the higher index to make room for the lost operons
                if index2 > index1:
                    targetIndex = index2
                else:
                    targetIndex = index1

                if targetIndex in arrangedFragments:
                    arrangedFragments[targetIndex].append(fragment)
                else:
                    arrangedFragments[targetIndex] = []
                    arrangedFragments[targetIndex].append(fragment)

    #Handle the remaining regions
    arrangedFragments = insertRegionIntoDictionary(TR, NFCR, arrangedFragments)
    arrangedFragments = insertRegionIntoDictionary(IR, NFCR, arrangedFragments)
    arrangedFragments = insertRegionIntoDictionary(ITR, NFCR, arrangedFragments)

    geneIndex = 0
    fragmentIndex = 0
    fragmentIndexes = arrangedFragments.keys()
    fragmentIndexes.sort()

    #Create the ancestral genome
    for x in range(0, len(fragmentIndexes)):
        key = fragmentIndexes[x]
        fragments = arrangedFragments[key]

        for fragment in fragments:
            originalSequence = ''
            negativeOrientation = False
            geneSequence = copy.deepcopy(fragment.ancestralOperonGeneSequence)

            if len(geneSequence) == 1 and fragment.fragmentDetails1.originalSequence != '< o >' and fragment.fragmentDetails1.originalSequence != '< t >':
                description = 'Singleton'
            else:
                description = 'Operon'

            if fragment.fragmentDetails1.isNegativeOrientation:
                negativeOrientation = True
                originalSequence = '-'

            originalSequence += '['
            for y in range(0, len(geneSequence)):
                gene = geneSequence[y]

                if y != (len(geneSequence) - 1):
                    originalSequence += gene + ', '
                else:
                    originalSequence += gene
            originalSequence += ']'

            newFragment = GenomeFragment(fragmentIndex, originalSequence, geneSequence, geneIndex, description, negativeOrientation)
            ancestralFragments.append(newFragment)

            index+=1
            geneIndex += len(geneSequence)

    return ancestralFragments

######################################################
# insertRegionIntoDictionary
# Parameters:
# Description: inserts region into fragment dictionary
######################################################
def insertRegionIntoDictionary(regions, NFCR, arrangedFragments):

    #Transposed/Inverted/Inverted Transposed regions
    for region in regions:
        count = 0

        #Iterate through all of the fragments in the region and count the number of times they appear in the neighbor's forward conserved region
        for x in range(0, len(region)):
            fragment = region[x]
            match = checkForFragment(fragment, NFCR)
            if match:
                count += 1
        #Insert the fragments into the dictionary based on the counter
        for x in range(0, len(region)):
            fragment = region[x]
            if count > 0:
                targetIndex = fragment.fragmentDetails1.fragmentIndex #Same arrangement exists in the neighbor
            else:
                targetIndex = fragment.fragmentDetails2.fragmentIndex #Neighbor's arrangement does not match Strain 1

            if targetIndex in arrangedFragments:
                arrangedFragments[targetIndex].append(fragment)
            else:
                arrangedFragments[targetIndex] = []
                arrangedFragments[targetIndex].append(fragment)

    return arrangedFragments

######################################################
# checkForFragment
# Parameters:
# Description: returns true if a matching fragment is found
######################################################
def checkForFragment(fragment, NFCR):
    match = False
    for region in NFCR:
        for x in range(0, len(region)):
            currFragment = region[x]
            if fragment.genome1Name == currFragment.genome1Name and fragment.fragmentDetails1.fragmentIndex == currFragment.fragmentDetails1.fragmentIndex:
                return True
    return match

######################################################
# determineAncestralFragmentArrangementWithoutNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the siblings
######################################################
def determineAncestralFragmentArrangementWithoutNeighbor(FCR, TR, IR, ITR, LR):
    ancestralFragments = []
    arrangedFragments = {}

    #Lost regions
    for x in range(0, len(LR)):
        fragment = LR[x]
        index = fragment.fragmentDetails1.fragmentIndex

        if index in arrangedFragments:
            arrangedFragments[index].append(fragment)
        else:
            arrangedFragments[index] = []
            arrangedFragments[index].append(fragment)

    #Forward conserved regions
    for region in FCR:
        for x in range(0, len(region)):
            fragment = region[x]
            index1 = fragment.fragmentDetails1.fragmentIndex
            index2 = fragment.fragmentDetails2.fragmentIndex

            if (index1 == index2): #Matching indexes
                if index1 in arrangedFragments:
                    arrangedFragments[index1].append(fragment)
                else:
                    arrangedFragments[index1] = []
                    arrangedFragments[index1].append(fragment)
            else: #Indexes don't match means there's a gap somewhere, insert into the higher index
                if index2 > index1:
                    targetIndex = index2
                else:
                    targetIndex = index1

                if targetIndex in arrangedFragments:
                    arrangedFragments[targetIndex].append(fragment)
                else:
                    arrangedFragments[targetIndex] = []
                    arrangedFragments[targetIndex].append(fragment)

    arrangedFragments = insertFragmentsIntoGenome(TR, arrangedFragments)
    arrangedFragments = insertFragmentsIntoGenome(IR, arrangedFragments)
    arrangedFragments = insertFragmentsIntoGenome(ITR, arrangedFragments)

    keyList = arrangedFragments.keys()
    keyList.sort()

    index = 0
    geneIndex = 0
    for x in range(0, len(keyList)):
        key = keyList[x]
        fragments = arrangedFragments[key]
        for fragment in fragments:
            originalSequence = ''
            negativeOrientation = False
            seq = copy.deepcopy(fragment.ancestralOperonGeneSequence)

            if len(seq) == 1 and fragment.fragmentDetails1.originalSequence != '< o >' and fragment.fragmentDetails1.originalSequence != '< t >':
                description = 'Singleton'
            else:
                description = 'Operon'

            if fragment.fragmentDetails1.isNegativeOrientation:
                negativeOrientation = True
                originalSequence = '-'

            originalSequence += '['
            for y in range(0, len(seq)):
                gene = seq[y]
                if y != (len(seq) - 1):
                    originalSequence += gene + ', '
                else:
                    originalSequence += gene
            originalSequence += ']'

            newFragment = GenomeFragment(index, originalSequence, seq, geneIndex, description, negativeOrientation)
            ancestralFragments.append(newFragment)

            index+=1
            geneIndex += len(seq)

    return ancestralFragments

######################################################
# insertFragmentsIntoGenome
# Parameters:
# Description: Inserts given fragments into genome
######################################################
def insertFragmentsIntoGenome(fragments, arrangedFragments):
    #Transposed/Inverted/Transposed Inverted
    for region in fragments:
        for x in range(0, len(region)):
            fragment = region[x]
            index1 = fragment.fragmentDetails1.fragmentIndex

            if index1 in arrangedFragments:
                arrangedFragments[index1].append(fragment)
            else:
                arrangedFragments[index1] = []
                arrangedFragments[index1].append(fragment)