import math
import copy
import globals
from GenomeFragment import GenomeFragment

###################################################
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
        
        index = 0
        
        oppositeOrientationCount = 0 #Counts the number of operons with opposite orientation within a consecutive region
        
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
        
        #Counts the number of mapped operons with the opposite orientation            
        if currFragment.fragmentDetails1.isNegativeOrientation != currFragment.fragmentDetails2.isNegativeOrientation:
            oppositeOrientationCount += 1

        while index < len(fragmentsCopy) and foundNeighbor and len(fragmentsCopy) > 0: #Continue looping as long as we keep on finding a neighbor
            foundNeighbor = False #Reset the tracker

            prevFragment = consecutiveRegion[len(consecutiveRegion)-1] #Gets the last operon in the the consecutive operons list
            currFragment = fragmentsCopy[index] #Get the next available operon
            yDistance = abs(currFragment.fragmentDetails2.point - prevFragment.fragmentDetails2.point) #The distance on the y-axis
            xDistance = abs(currFragment.fragmentDetails1.point - prevFragment.fragmentDetails1.point) #The distance on the x-axis
            
            #cMainDiagonalDistance = abs(currFragment.fragmentDetails1.point - currFragment.fragmentDetails2.point) #Current points distance from the main diagonal
            #pMainDiagonalDistance = abs(prevFragment.fragmentDetails1.point - prevFragment.fragmentDetails2.point) #Previous points distance from the main diagonal
            #not(pMainDiagonalDistance == 0 and cMainDiagonalDistance > 0)
            if yDistance < globals.yDistanceThreshold and xDistance < globals.xDistanceThreshold: #If the y-Distance is less than the threshold add it to the consecutive region
                consecutiveRegion.append(fragmentsCopy.pop(index))
                foundNeighbor = True #Indicates we found a consecutive region
                #skipCounter = 0 #Reset the skip counter
                
                if currFragment.fragmentDetails1.isNegativeOrientation != currFragment.fragmentDetails2.isNegativeOrientation:
                    oppositeOrientationCount += 1
                
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
        sameOrientation = len(consecutiveRegion) - oppositeOrientationCount # number of positive oriented operons
        if oppositeOrientationCount > 0 and oppositeOrientationCount >= sameOrientation:
            #This is either an inversion or an inverted transposition
            if aboveMainDiagonal == True and belowMainDiagonal == True and yIncreaseCounter < yDecreaseCounter:
                invertedRegions.append(consecutiveRegion) #Crosses the main diagonal and y is decreasing
            else:
                invertedTransposedRegions.append(consecutiveRegion) #If does not cross main diagonal and y is not decreasing, treat as inverted transposition
        else:
            #This is either a forward conserved region or a transposition
            if (minMainDiagonalDistance < 3 and yIncreaseCounter > yDecreaseCounter) or (minMainDiagonalDistance == 0 and len(consecutiveRegion) == 1): #If the distance from the main diagonal is 0 and y is increasing
                conservedForwardRegions.append(consecutiveRegion)
            else:
                transposedForwardRegions.append(consecutiveRegion)
                
    if globals.printToConsole:
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
                print('%s, %s : %s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point, region[x].fragmentDetails1.startPositionInGenome, region[x].fragmentDetails2.startPositionInGenome))
        for region in transposedForwardRegions:
            print('Forward Transposed Region')
            for x in range(0, len(region)):
                print('%s, %s : %s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point, region[x].fragmentDetails1.startPositionInGenome, region[x].fragmentDetails2.startPositionInGenome))
        for region in invertedRegions:
            print('Inverted Region')
            for x in range(0, len(region)):
                print('%s, %s : %s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point, region[x].fragmentDetails1.startPositionInGenome, region[x].fragmentDetails2.startPositionInGenome))
        for region in invertedTransposedRegions:
            print('Inverted Transposed Region')
            for x in range(0, len(region)):
                print('%s, %s : %s, %s' %(region[x].fragmentDetails1.point, region[x].fragmentDetails2.point, region[x].fragmentDetails1.startPositionInGenome, region[x].fragmentDetails2.startPositionInGenome))

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
    
    if globals.printToConsole:
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
            temp1 = temp1[0:(len(temp1) - 2)]
            temp1 += ';'

            for y in range(0, len(seq2)):
                currPos = startPos2 + y
                currGene = seq2[y]
                temp2 += currGene + ' ' + str(currPos) + ', '
            temp2 = temp2[0:(len(temp2) - 2)]
            temp2 += ';'
        temp1 = temp1.strip() + '|' #End of region
        temp2 = temp2.strip() + '|' #End of region

    return temp1, temp2

######################################################
# determineAncestralFragmentArrangementUsingNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the siblings and the neighbor
######################################################
def determineAncestralFragmentArrangementUsingNeighbor(FCR, TR, IR, ITR, LR, NFCR, NTR, NIR, NITR, NLR, strain1, strain2):
    #Initialize
    arrangedFragments = {}

    #Insert the lost operons at the positions they were lost in
    for fragment in LR:
        index = fragment.fragmentDetails1.fragmentIndex
        fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon

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
                fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                if index1 in arrangedFragments:
                    arrangedFragments[index1].append(fragment)
                else:
                    arrangedFragments[index1] = []
                    arrangedFragments[index1].append(fragment)
            else: #The positions have not been conserved, means there's some operons lost between the two genomes. Insert into the higher index to make room for the lost operons
                if index2 > index1:
                    targetIndex = index2
                    fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails2.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                else:
                    targetIndex = index1
                    fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon

                if targetIndex in arrangedFragments:
                    arrangedFragments[targetIndex].append(fragment)
                else:
                    arrangedFragments[targetIndex] = []
                    arrangedFragments[targetIndex].append(fragment)

    #Transposed regions
    arrangedFragments, details1, details1Counter, details2, details2Counter = insertRegionIntoDictionary(TR, NFCR, arrangedFragments)
    if len(details1Counter) > 0:
        for size, count in details1Counter.items():
            if size in strain2.transpositionCounts:
                strain2.transpositionCounts[size] += count
            else:
                strain2.transpositionCounts[size] = count
        strain2.transpositionDetails += details1
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain1.transpositionCounts:
                strain1.transpositionCounts[size] += count
            else:
                strain1.transpositionCounts[size] = count
        strain1.transpositionDetails += details2

    #Inverted regions
    arrangedFragments, details1, details1Counter, details2, details2Counter = insertRegionIntoDictionary(IR, NFCR, arrangedFragments)
    if len(details1Counter) > 0:
        for size, count in details1Counter.items():
            if size in strain2.inversionCounts:
                strain2.inversionCounts[size] += count
            else:
                strain2.inversionCounts[size] = count
        strain2.inversionDetails += details1
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain1.inversionCounts:
                strain1.inversionCounts[size] += count
            else:
                strain1.inversionCounts[size] = count
        strain1.inversionDetails += details2

    #Inverted Transposed regions
    arrangedFragments, details1, details1Counter, details2, details2Counter = insertRegionIntoDictionary(ITR, NFCR, arrangedFragments)
    if len(details1Counter) > 0:
        for size, count in details1Counter.items():
            if size in strain2.invertedTranspositionCounts:
                strain2.invertedTranspositionCounts[size] += count
            else:
                strain2.invertedTranspositionCounts[size] = count
        strain2.invertedTranspositionDetails += details1
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain1.invertedTranspositionCounts:
                strain1.invertedTranspositionCounts[size] += count
            else:
                strain1.invertedTranspositionCounts[size] = count
        strain1.invertedTranspositionDetails += details2

    #Construct the return the genome
    ancestralFragments = constructGenome(arrangedFragments)
    return ancestralFragments, strain1, strain2

######################################################
# insertRegionIntoDictionary
# Parameters:
# Description: inserts region into fragment dictionary
######################################################
def insertRegionIntoDictionary(regions, NFCR, arrangedFragments):
    details1 = ''
    details1Counter = {}
    details2 = ''
    details2Counter = {}

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
        addedDetails1 = False
        addedDetails2 = False
        size1 = 0
        size2 = 0
        
        #Sort region so the correct sequence is generated
        if count > 0:
            region.sort(key=lambda x:x.fragmentDetails2.fragmentIndex, reverse=False)
        else:
            region.sort(key=lambda x:x.fragmentDetails1.fragmentIndex, reverse=False)
        
        for x in range(0, len(region)):
            fragment = region[x]
            if count > 0:
                targetIndex = fragment.fragmentDetails1.fragmentIndex #Neighbor's arrangement does not match Strain 1
                fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                #Constructs the description of the region
                seq = fragment.fragmentDetails2.sequence
                size1 += len(seq)
                startPosition = fragment.fragmentDetails2.startPositionInGenome
                temp = ''
                for i in range(0, len(seq)):
                    if fragment.fragmentDetails2.isNegativeOrientation:
                        temp = seq[i] + ' ' + str(startPosition + len(seq) - i - 1) + ', ' + temp
                    else:
                        temp += seq[i] + ' ' + str(startPosition + i) + ', ' 
                temp = temp[0:(len(temp) - 2)]
                temp += '; '
                details1 += temp.strip()
                addedDetails1 = True
                
            else:
                targetIndex = fragment.fragmentDetails2.fragmentIndex #Same arrangement exists in the neighbor
                fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails2.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                #Constructs the description of the region
                seq = fragment.fragmentDetails1.sequence
                size2 += len(seq)
                startPosition = fragment.fragmentDetails1.startPositionInGenome
                temp = ''
                for i in range(0, len(seq)):
                    if fragment.fragmentDetails1.isNegativeOrientation:
                        temp = seq[i] + ' ' + str(startPosition + len(seq) - i - 1) + ', ' + temp
                    else:
                        temp += seq[i] + ' ' + str(startPosition + i) + ', '
                temp = temp[0:(len(temp) - 2)]
                temp += '; '
                details2 += temp.strip()
                addedDetails2 = True

            if targetIndex in arrangedFragments:
                arrangedFragments[targetIndex].append(fragment)
            else:
                arrangedFragments[targetIndex] = []
                arrangedFragments[targetIndex].append(fragment)

        #Add a delimiter if a region was added and increment the appropriate counter
        if addedDetails2:
            details2 = details2.strip() + '|'
            if size2 in details2Counter:
                details2Counter[size2] += 1
            else:
                details2Counter[size2] = 1
        if addedDetails1:
            details1 = details1.strip() + '|'
            if size1 in details1Counter:
                details1Counter[size1] += 1
            else:
                details1Counter[size1] = 1

    return arrangedFragments, details1, details1Counter, details2, details2Counter

######################################################
# checkForFragment
# Parameters:
# Description: returns true if a matching fragment is found
######################################################
def checkForFragment(fragment, NFCR):
    match = False
    for region in NFCR:
        filteredList = iter(filter(lambda x:x.fragmentDetails1.fragmentIndex == fragment.fragmentDetails1.fragmentIndex, region))
        foundFragment = next(filteredList, None)
        if foundFragment:
            return True
    return match

######################################################
# determineAncestralFragmentArrangementWithoutNeighbor
# Parameters:
# Description: Orders the ancestral fragments based on the siblings
######################################################
def determineAncestralFragmentArrangementWithoutNeighbor(FCR, TR, IR, ITR, LR, strain):
    arrangedFragments = {}

    #Insert the lost operons at the positions they were lost in
    for fragment in LR:
        index = fragment.fragmentDetails1.fragmentIndex
        fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon
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
                fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                if index1 in arrangedFragments:
                    arrangedFragments[index1].append(fragment)
                else:
                    arrangedFragments[index1] = []
                    arrangedFragments[index1].append(fragment)
            else: #The positions have not been conserved, means there's some operons lost between the two genomes. Insert into the higher index to make room for the lost operons
                if index2 > index1:
                    targetIndex = index2
                    fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails2.isNegativeOrientation) #Identifies the orientation of the ancestral operon
                else:
                    targetIndex = index1
                    fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon

                if targetIndex in arrangedFragments:
                    arrangedFragments[targetIndex].append(fragment)
                else:
                    arrangedFragments[targetIndex] = []
                    arrangedFragments[targetIndex].append(fragment)
    #Transpositions
    arrangedFragments, details2, details2Counter = insertFragmentsIntoGenome(TR, arrangedFragments)
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain.transpositionCounts:
                strain.transpositionCounts[size] += count
            else:
                strain.transpositionCounts[size] = count
        strain.transpositionDetails += details2

    #Inversions
    arrangedFragments, details2, details2Counter = insertFragmentsIntoGenome(IR, arrangedFragments)
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain.inversionCounts:
                strain.inversionCounts[size] += count
            else:
                strain.inversionCounts[size] = count
        strain.inversionDetails += details2

    #Inverted Transpositions
    arrangedFragments, details2, details2Counter = insertFragmentsIntoGenome(ITR, arrangedFragments)
    if len(details2Counter) > 0:
        for size, count in details2Counter.items():
            if size in strain.invertedTranspositionCounts:
                strain.invertedTranspositionCounts[size] += count
            else:
                strain.invertedTranspositionCounts[size] = count
        strain.invertedTranspositionDetails += details2

    #Construct and return the genome
    ancestralFragments = constructGenome(arrangedFragments)
    return ancestralFragments, strain

######################################################
# addDetails
# Parameters:
# Description: Appends inversion, transposition, inverted transposition details to a dictionary and String
######################################################
def addDetails(allDetails, allDetailsCounter, detailsCounter, details):
    if len(detailsCounter) > 0:
        for size, count in detailsCounter.items():
            if size in allDetailsCounter:
                allDetailsCounter[size] += count
            else:
                 allDetailsCounter[size] = 1
        allDetails += details
    return allDetails, allDetailsCounter

######################################################
# insertFragmentsIntoGenome
# Parameters:
# Description: Inserts given fragments into genome
######################################################
def insertFragmentsIntoGenome(fragments, arrangedFragments):
    details2 = ''
    details2Counter = {}

    #Transposed/Inverted/Transposed Inverted
    for region in fragments:
        size = 0
        addedDetails2 = False #Tracks whether we added a region
        
        #Sort region
        region.sort(key=lambda x:x.fragmentDetails2.fragmentIndex, reverse=False)
        for x in range(0, len(region)):
            fragment = region[x]
            index1 = fragment.fragmentDetails1.fragmentIndex
            fragment.setAncestralOperonNegativeOrientation(fragment.fragmentDetails1.isNegativeOrientation) #Identifies the orientation of the ancestral operon

            #We will assume Strain 2 was transposed, inverted, inverted transposed if no neighbor is present
            seq = fragment.fragmentDetails2.sequence
            size += len(seq)
            startPosition = fragment.fragmentDetails2.startPositionInGenome
            temp = ''
            for i in range(0, len(seq)):
                
                if fragment.fragmentDetails2.isNegativeOrientation:
                    temp = seq[i] + ' ' + str(startPosition + len(seq) - i - 1) + ', ' + temp
                else:
                    temp += seq[i] + ' ' + str(startPosition + i) + ', ' 

            temp = temp[0:(len(temp) - 2)]
            temp += '; '
            details2 += temp.strip()
            addedDetails2 = True

            if index1 in arrangedFragments:
                arrangedFragments[index1].append(fragment)
            else:
                arrangedFragments[index1] = []
                arrangedFragments[index1].append(fragment)
        #Add a delimiter if a region was added and increment the appropriate counter
        if addedDetails2:
            details2 = details2.strip() + '|'
            if size in details2Counter:
                details2Counter[size] += 1
            else:
                details2Counter[size] = 1

    return arrangedFragments, details2, details2Counter

######################################################
# constructGenome
# Parameters:
# Description: Constructs the ancestral genome
######################################################
def constructGenome(arrangedFragments):
    #Initialize variables
    geneIndex = 0
    fragmentIndex = 0
    ancestralFragments = []
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
            positiveOrientationGeneSequence = copy.deepcopy(fragment.ancestralOperonGeneSequence)

            if len(geneSequence) == 1 and fragment.fragmentDetails1.sequence[0] != '< o >' and fragment.fragmentDetails1.sequence[0] != '< t >':
                description = 'Singleton'
            else:
                #Make sure description is correct
                if fragment.fragmentDetails1.sequence[0] == '< o >':
                    description = 'Origin'
                elif fragment.fragmentDetails1.sequence[0] == '< t >':
                    description = 'Terminus'
                else:
                    description = 'Operon'

            if fragment.ancestralOperonNegativeOrientation == True:
                negativeOrientation = True
                originalSequence = '-'
                geneSequence.reverse()
            
            #Brackets only on the oerons with length > 1
            if len(geneSequence) != 1:
                originalSequence += '['
            
            for y in range(0, len(geneSequence)):
                gene = geneSequence[y]

                if y != (len(geneSequence) - 1):
                    originalSequence += gene + ', '
                else:
                    originalSequence += gene
            #Brackets only on operons with length > 1
            if len(geneSequence) != 1:
                originalSequence += ']'

            newFragment = GenomeFragment(fragmentIndex, originalSequence, positiveOrientationGeneSequence, geneIndex, description, negativeOrientation)
            ancestralFragments.append(newFragment)
            
            #Update the deletion details if needed
            if len(fragment.deletionDetailsList) > 0:
                for item in fragment.deletionDetailsList:
                    item.ancestralFragmentId = fragmentIndex
                    newFragment.deletionDetailsList.append(item)

            fragmentIndex += 1
            geneIndex += len(geneSequence)

    #Make sure the fragments are sorted by the index
    ancestralFragments.sort(key=lambda x:x.fragmentIndex, reverse=False)
    return ancestralFragments