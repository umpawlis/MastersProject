import time
import os.path
from Bio import Phylo
from strain import Strain
from GlobalAlignmentModule import findOrthologsByGlobalAlignment
from LocalAlignmentModule import findOrthologsByLocalAlignment
from SelfGlobalAlignmentModule import findOrthologsBySelfGlobalAlignment
from FileService import createFile
from FileService import appendToFile
import globals
import matplotlib.pyplot as plt
import numpy as np
import copy

#Parameters that user will pass in
newickFileName = 'Bacillus_Tree.dnd'

#Global variables used in script
strains = [] #Global variable that stores all of the strains in the phylogeny
ancestralCounter = 0 #Counter used to create a unique name for the ancestor
deletionCost = 1
substitutionCost = 1
codonCost = 0.5
outputFile = 'ApplicationOutput.txt'

#################################################
########Functions used in this script############
#################################################

######################################################
# createBarGraph
# Parameters:
# Description:
######################################################
def createBarGraph(dictionary, title):

    if dictionary != None and len(dictionary) > 0:
        keys = list(dictionary.keys())
        keys.sort()

        y_pos = np.arange(len(keys))

        performance = []
        for key in keys:
            performance.append(dictionary[key])

        plt.bar(y_pos, performance, align='center', alpha=0.5)
        plt.xticks(y_pos, keys)
        plt.ylabel('Number of Occurrences')
        plt.xlabel('Size of Occurrence')
        plt.title(title)
        plt.show()

######################################################
# createDotPlot
# Parameters:
# Description:
######################################################
def createDotPlot(events, strain1, strain2):
    #Stores all of the coordinates
    x_coord = []
    y_coord = []

    #The green ones represent operons with no differences
    green_x_coord = []
    green_y_coord = []

    #Yellow ones represent scores between 1 and 2
    yellow_x_coord = []
    yellow_y_coord = []

    #Red ones represent scores between 3 and above
    orange_x_coord = []
    orange_y_coord = []

    #Blue ones represent a local alignment
    red_x_coord = []
    red_y_coord = []

    print("x" * 70)
    for i in range(0, len(events)):
        if events[i].technique == 'Global Alignment' or events[i].technique == 'Local Alignment':
            #Assign the coords to the appropriate array
            if events[i].technique == 'Local Alignment':
                red_x_coord.append(events[i].operon1Index)
                red_y_coord.append(events[i].operon2Index)
            elif events[i].score == 0:
                green_x_coord.append(events[i].operon1Index)
                green_y_coord.append(events[i].operon2Index)
            elif events[i].score == 1 or events[i].score == 2:
                yellow_x_coord.append(events[i].operon1Index)
                yellow_y_coord.append(events[i].operon2Index)
            else:
                orange_x_coord.append(events[i].operon1Index)
                orange_y_coord.append(events[i].operon2Index)

            #Get all coordinates into a single array
            x_coord.append(events[i].operon1Index)
            y_coord.append(events[i].operon2Index)

            #print('x-axis: %s, y-axis: %s' %(trackingEvents[i].getGenome1OperonIndex(), trackingEvents[i].getGenome2OperonIndex()))

    #If we have any coordinates to plot, display them
    if len(green_x_coord) > 0 or len(yellow_x_coord) > 0 or len(orange_x_coord) > 0 or len(red_x_coord) > 0:
        f = plt.figure()
        plt.title("Orthologous Operon Mapping")
        plt.plot(green_x_coord, green_y_coord, 'o', color = 'green')
        plt.plot( yellow_x_coord, yellow_y_coord, 'o', color = 'yellow')
        plt.plot(orange_x_coord, orange_y_coord, 'o', color = 'orange')
        plt.plot(red_x_coord, red_y_coord, 'o', color = 'red')
        plt.axis([0, len(events)+5, 0, len(events)+5])
        plt.ylabel('Operon Position in %s' % (strain1.getName()))
        plt.xlabel('Operon Position in %s' % (strain2.getName()))
        plt.show()
        f.savefig("%s %s.pdf" %(strain1.getName(), strain2.getName()), bbox_inches='tight')
    else:
        print('No plot to display!')
    print("x" * 70)

######################################################
# formatAndOrientOperon
# Parameters:
# Description: Formats the operon into a string with the correct orientation
######################################################
def formatAndOrientOperon(events, neighborEvents):
    if events != None and len(events) > 0:
        for event in events:
            if event.reversedOp1 == False and event.reversedOp2 == False: #Operons can be left as is
                print('These operons were not reversed')
                if event.originallyNegativeOrientationOp1:
                    stringAncestralOperon = formatOperonIntoString(event, True)
                    event.setStringAncestralOperon(stringAncestralOperon)
                    print('This is the - ancestral operon reconstucted %s' %(stringAncestralOperon))
                else:
                    stringAncestralOperon = formatOperonIntoString(event, False)
                    event.setStringAncestralOperon(stringAncestralOperon)
                    print('This is the + ancestral operon reconstucted %s' %(stringAncestralOperon))
            else:
                print('These operons have opposite orientation!') #Need to check neighbor to determine orientation
                filteredList = iter(filter(lambda x: x.genome1Operon == event.genome1Operon, neighborEvents))
                neighborEvent = next(filteredList, None)
                if neighborEvent != None:
                    print('Found same operon in neighbor')
                    if neighborEvent.originallyNegativeOrientationOp2 == True and (event.originallyNegativeOrientationOp1 == True or event.originallyNegativeOrientationOp2 == True):
                        stringAncestralOperon = formatOperonIntoString(event, True) #Neighbor operon and 1 other operon is - so - majority wins
                        event.setStringAncestralOperon(stringAncestralOperon)
                        print('This is the - ancestral operon reconstucted %s' %(stringAncestralOperon))
                    else:
                        stringAncestralOperon = formatOperonIntoString(event, False) #Neighbor operon and 1 other operon is + so + majority wins
                        event.setStringAncestralOperon(stringAncestralOperon)
                        print('This is the + ancestral operon reconstucted %s' %(stringAncestralOperon))
                else:
                    print('Operon not found in neighbor strain') #No matching operon found in neighbor, just assign an orientation
                    stringAncestralOperon = formatOperonIntoString(event, True)
                    event.setStringAncestralOperon(stringAncestralOperon)
                    print('This is the - ancestral operon reconstucted %s' %(stringAncestralOperon))
    return events

######################################################
# formatOperonIntoString
# Parameters:
# Description: Formats the operon into either a postive or negative oriented string
######################################################
def formatOperonIntoString(event, isNegativeOrientation):
    ancestralOperonString = ''
    genes = event.ancestralOperonGeneSequence

    if len(genes) == 1: #We are dealing with a singleton
        if isNegativeOrientation:
            print('This is a singleton in the - orientation')
            ancestralOperonString = '-' + str(genes[0])
        else:
            print('This is a singleton in the + orientation')
            ancestralOperonString = str(genes[0])
    else: #We are dealing with an operon
        if isNegativeOrientation:
            print('This is an operon in the - orientation')
            ancestralOperonString = '-['
            for x in range(0, len(genes)):
                if x == len(genes) - 1:
                    ancestralOperonString += genes[x] + ']'
                else:
                    ancestralOperonString += genes[x] + ', '
        else:
            print('This is an operon in the + orientation')
            ancestralOperonString = '['
            for x in range(0, len(genes)):
                if x == len(genes) - 1:
                    ancestralOperonString += genes[x] + ']'
                else:
                    ancestralOperonString += genes[x] + ', '
    return ancestralOperonString

######################################################
# computeOperonArrangements
# Parameters:
# Description: Takes a list of events and computes Forward Conserved, Transposed Forward, Inverted, Inverted Transposed, and Lost regions
######################################################
def computeOperonArrangements(events):
    conservedForwardRegions = []
    transposedForwardRegions = []
    invertedRegions = []
    invertedTransposedRegions = []
    lostRegions = []
    
    #Sort events by x coordinates
    eventsCopy = copy.deepcopy(events)
    eventsCopy.sort(key=lambda x : x.operon1Index, reverse=False)
    
    while len(eventsCopy) > 0:
        currEvent = eventsCopy.pop(0)
        
        if currEvent.operon1Index == -1 or currEvent.operon2Index == -1: #We are dealing with a lost operon
            lostRegions.append(currEvent)
        else:
            consecutiveRegion = []
            consecutiveRegion.append(currEvent)
            foundNeighbor = True
            yIncreaseCount = 0
            yDecreaseCount = 0
            aboveMainDiagonal = False
            belowMainDiagonal = False
            minMainDiagonalDistance = abs(currEvent.operon1Index - currEvent.operon2Index)
            if currEvent.operon1Index - currEvent.operon2Index < 0:
                aboveMainDiagonal = True
            if currEvent.operon1Index - currEvent.operon2Index > 0:
                belowMainDiagonal = True
            while foundNeighbor:
                foundNeighbor = False
                if len(eventsCopy) > 0:
                    prevEvent = consecutiveRegion[len(consecutiveRegion)-1]
                    currEvent = eventsCopy[0]
                    distance = abs(currEvent.operon2Index - prevEvent.operon2Index)
                    if distance < globals.yDistanceThreshold:
                        consecutiveRegion.append(eventsCopy.pop(0))
                        foundNeighbor = True
                        currMainDiagonalDistance = abs(currEvent.operon1Index - currEvent.operon2Index)
                        #Tracks whether this is a transposition
                        if currMainDiagonalDistance < minMainDiagonalDistance:
                            minMainDiagonalDistance = currMainDiagonalDistance
                        #Tracks whether we cross the main diagonal
                        if currEvent.operon1Index - currEvent.operon2Index < 0:
                            aboveMainDiagonal = True
                        if currEvent.operon1Index - currEvent.operon2Index > 0:
                            aboveMainDiagonal = True
                        #Tracks if the points are moving up or down
                        if prevEvent.operon2Index < currEvent.operon2Index:
                            yIncreaseCount += 1
                        else:
                            yDecreaseCount += 1
            if yIncreaseCount > yDecreaseCount or len(consecutiveRegion) == 1:
                if minMainDiagonalDistance < globals.yDistanceThreshold:
                    conservedForwardRegions.append(consecutiveRegion)
                else:
                    transposedForwardRegions.append(consecutiveRegion)
            else:
                if aboveMainDiagonal == True and belowMainDiagonal == True:
                    invertedRegions.append(consecutiveRegion)
                else:
                    invertedTransposedRegions.append(consecutiveRegion)
    print('Stats:')
    print('Total number of tracking events: %s' % (len(events)))
    print('Total number of forward conserved regions: %s' % len(conservedForwardRegions))
    print('Total number of forward transposed regions: %s' % len(transposedForwardRegions))
    print('Total number of inverted regions: %s' % len(invertedRegions))
    print('Total number of inverted transposed regions: %s' % len(invertedTransposedRegions))
    print('Total number of lost operons: %s' % len(lostRegions))

    for region in conservedForwardRegions:
        print('Forward Conserved Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].operon1Index, region[x].operon2Index))
    for region in transposedForwardRegions:
        print('Forward Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].operon1Index, region[x].operon2Index))
    for region in invertedRegions:
        print('Inverted Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].operon1Index, region[x].operon2Index))
    for region in invertedTransposedRegions:
        print('Inverted Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].operon1Index, region[x].operon2Index))
    for operon in lostRegions:
        print('Lost operon')
        print('%s, %s' %(operon.operon1Index, operon.operon2Index))
        
    return conservedForwardRegions, transposedForwardRegions, invertedRegions, invertedTransposedRegions, lostRegions

######################################################
# processStrains
# Parameters:
# Description: Takes two related strains and a close neighbor and constructs the events for both comparisons
######################################################
def processStrains(strain1, strain2, neighborStrain):
    ancestralSequence = []
    events = []

    print('Constructing events of the following siblings: %s, %s' %(strain1.getName(), strain2.getName()))
    events = constructEvents(strain1, strain2)
    temp1 = copy.deepcopy(globals.sizeDuplications)
    temp2 = copy.deepcopy(globals.sizeDeletions)
    
    print('Constructing dot plot for the following siblings: %s, %s' %(strain1.getName(), strain2.getName()))
    createDotPlot(events, strain1, strain2)
    createBarGraph(strain1.duplicationSizes, 'Distribution of Duplications for %s'%(strain1.getName()))
    createBarGraph(strain2.duplicationSizes, 'Distribution of Duplications for %s'%(strain2.getName()))
    createBarGraph(strain1.deletionSizes, 'Distribution of Deletions for %s'%(strain1.getName()))
    createBarGraph(strain2.deletionSizes, 'Distribution of Deletions for %s'%(strain2.getName()))
    #createBarGraph(events, strain1, strain2, globals.localSizeDuplications, 'Distribution of Duplications %s vs %s' % (strain1.getName(), strain2.getName()))
    #createBarGraph(events, strain1, strain2, globals.localSizeDeletions, 'Distribution of Deletions %s vs %s' % (strain1.getName(), strain2.getName()))
    
    neighborEvents = []
    if neighborStrain != None:#Had to put this check in since it was causing it crash for root
        neighborEvents = constructEvents(copy.deepcopy(strain1), neighborStrain)
    events = formatAndOrientOperon(events, neighborEvents) #Inserts into each event an ancestral operon in a string format

    FCR, TFCR, IR, ITR, LR = computeOperonArrangements(events)
    regionCount = len(TFCR) + len(IR) + len(ITR)
    if regionCount != 0 and len(neighborEvents) != 0:
        NFCR, NTFCR, NIR, NITR, NLR = computeOperonArrangements(neighborEvents)
        ancestralSequence = constructGenome(FCR, TFCR, IR, ITR, LR, NFCR, NTFCR, NIR, NITR, NLR)
    else: #Sort events by x coordinate and reconstruct the genome
        if regionCount == 0:
            print('No inverted or transposed regions!')
        elif len(neighborEvents):
            print('No neighboring events!')
        eventsCopy = copy.deepcopy(events)
        eventsCopy.sort(key=lambda x : x.operon1Index, reverse=False)
        for event in eventsCopy:
            ancestralSequence.append(event.stringAncestralOperon)

    #Clear out the global trackers since they'll have comparisons from the neighbor alignment
    globals.localSizeDuplications.clear()
    globals.localSizeDeletions.clear()
    globals.sizeDuplications.clear()
    globals.sizeDeletions.clear()

    #Set the global trackers of the phylogeny to what it was when we compared the siblings b/c at that point the dictionary is not messed up with neighbor data
    globals.sizeDuplications = temp1
    globals.sizeDeletions = temp2
    
    #Add content to output file
    appendStrainToFile(strain1)
    appendStrainToFile(strain2)
    
    return ancestralSequence, events

######################################################
# appendStrainToFile
# Parameters:
# Description: Adds content to the output file based on the strain's data
######################################################
def appendStrainToFile(strain):
    appendToFile(outputFile, "%s\n" % (strain.getName())) #Name of the strain
    
    temp = "Deletions;"
    if strain.deletionSizes != None and len(strain.deletionSizes):
        for size, count in strain.deletionSizes.items():
            temp += str(size) + ":" + str(count) + ","
    appendToFile(outputFile, "%s\n" %(temp)) #Deletions

######################################################
# constructGenome
# Parameters:
# Description: Takes regions detected from the siblings and compares and arranges those regions based on the neighbor
######################################################
def constructGenome(CFR, TFR, IR, ITR, LO, NCFR, NTFR, NIR, NITR, NLO):

    ancestralSequence = []
    dictionary = {}
    
    if LO != None and len(LO) > 0: #Deal with the lost operons
        for operon in LO:
            if operon.operon1Index in dictionary:
                dictionary[operon.operon1Index].append(operon)
            else:
                dictionary[operon.operon1Index] = [operon]
                
    if CFR != None and len(CFR) > 0: #Deal with the Forward Conserved Regions
        for region in CFR:
            for x in range(0, len(region)):
                if region[x].operon1Index in dictionary:
                    dictionary[region[x].operon1Index].append(region[x])
                else:
                    dictionary[region[x].operon1Index] = [region[x]]
    #Handle Transposed Forward Regions
    dictionary = processRegion(TFR, dictionary, NTFR, NIR, NITR)
    #Handle Inversed Regions
    dictionary = processRegion(IR, dictionary, NTFR, NIR, NITR)
    #Handle Inversed Transposed Regions
    dictionary = processRegion(ITR, dictionary, NTFR, NIR, NITR)
    #Now add the operons from the dictionary to the ancestral sequence
    keys = dictionary.keys()
    keys.sort()
    for key in keys:
        operons = dictionary[key]
        if not(operons is None) and len(operons) > 0:
            for operon in operons:
                ancestralSequence.append(operon.stringAncestralOperon)
    return ancestralSequence

######################################################
# processRegion
# Parameters: list of regions, dictionary which stores the position of the regions, neighbors regions
# Description: Takes an array of regions and checks if the postion of that region is conserved in the neighbor and adds them to the dictionary appropriately
######################################################
def processRegion(regions, dictionary, NTFR, NIR, NITR):
    if regions != None and len(regions) > 0:
        for region in regions:
            count = 0
            for x in range(0, len(region)):
                currEvent = region[x]
                #need to check neighbor if the same exists
                found = checkNeighorRegion(currEvent.operon1Index, NTFR)
                if found == True:
                    count += 1
                else:
                    found = checkNeighorRegion(currEvent.operon1Index, NIR)
                    if found == True:
                        count += 1
                    else:
                        found = checkNeighorRegion(currEvent.operon1Index, NITR)
                        if found:
                            count += 1
            if count > 0:
                for x in range (0, len(region)):
                    if region[x].operon1Index in dictionary:
                        dictionary[region[x].operon1Index].append(region[x])
                    else:
                        dictionary[region[x].operon1Index] = [region[x]]
            else:
                for x in range (0, len(region)):
                    if region[x].operon2Index in dictionary:
                        dictionary[region[x].operon2Index].append(region[x])
                    else:
                        dictionary[region[x].operon2Index] = [region[x]]
    return dictionary

######################################################
# checkNeighorRegion
# Parameters: index of operon, list of regions
# Description: Checks if the given operon(index) exists in a given array of regions
######################################################
def checkNeighorRegion(index, arrayOfRegions):
    if arrayOfRegions != None and len(arrayOfRegions) > 0:
        for region in arrayOfRegions:
            for x in range(0, len(region)):
                currEvent = region[x]
                if currEvent.operon1Index == index:
                    return True
    return False

######################################################
#countRemainingOperons
#Parameters: tracker - array of booleans that tracks whether an operon was paired (ortholog)
#Description: Takes an array of booleans and returns a count of the number of False values in the list
######################################################
def countRemainingOperons(tracker):
    count = 0
    for x in range(0, len(tracker)):
        if tracker[x] == False:
            count += 1
    return count

######################################################
# constructEvents
# Parameters:
# Description: Constructs the tracking events between two provided strains
######################################################
def constructEvents(strain1, strain2):
    globals.localSizeDuplications.clear()
    globals.localSizeDeletions.clear()

    events = []
    coverageTracker1 = {}
    coverageTracker2 = {}
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()

    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False

    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False

    #Global Alignment operation
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter = findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)
    print('Number of orthologous operons identified using Global Alignment %s' % (globalAlignmentCounter))

    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))

    #Local Alignment operation
    if numRemainingOperons1 > 0 and numRemainingOperons2 > 0:
        localAlignmentEvents, coverageTracker1, coverageTracker2, localAlignmentCounter = findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
        print('Number of orthologous operons identified using Local Alignment %s' % (localAlignmentCounter))

        numRemainingOperons1 = countRemainingOperons(coverageTracker1)
        numRemainingOperons2 = countRemainingOperons(coverageTracker2)
        print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))
        if len(localAlignmentEvents) > 0:
            events.extend(localAlignmentEvents)

    #Self Global Alignment
    if numRemainingOperons1 > 0:
        duplicationEvents1, lossEvents1, coverageTracker1 = findOrthologsBySelfGlobalAlignment(strain1, coverageTracker1)
        print('%s, duplicates identified %s and losses identified %s' % (strain1.getName(), len(duplicationEvents1), len(lossEvents1)))
        if len(lossEvents1) > 0:
            events.extend(lossEvents1)
    if numRemainingOperons2 > 0:
        duplicationEvents2, lossEvents2, coverageTracker2 = findOrthologsBySelfGlobalAlignment(strain2, coverageTracker2)
        print('%s, duplicates identified %s and losses identified %s' % (strain2.getName(), len(duplicationEvents2), len(lossEvents2)))
        if len(lossEvents2) > 0:
            events.extend(lossEvents2)

    #Verify there's no unmarked operons at this point
    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    if numRemainingOperons1 > 0 or numRemainingOperons2 > 0:
        print('Error! There are unmarked operons remaining!')

    return events

######################################################
# traverseNewickTree
# Parameters: node - Strain being currently processed, parentNode -  direct ancestor of node
# Description: Traverses a provided newick tree in post order traversal
######################################################
def traverseNewickTree(node, parentNode):
    global strains
    global ancestralCounter

    leftSibling = None
    rightSibling = None

    #Check there's a descendant
    if len(node.clades) > 0:
        leftSibling = traverseNewickTree(node.clades[0], node)
        if len(node.clades) > 1:
            rightSibling = traverseNewickTree(node.clades[1], node)

    #Retrieve strain or create one from the data file
    if not(node.name == None) and len(node.name) > 0:

        filteredList = iter(filter(lambda x: x.name == node.name, strains))
        foundStrain = next(filteredList, None)

        if (foundStrain != None):
            print('Retrieving strain from strains list: %s' % (foundStrain.getName()))
            return foundStrain
        else:
            newStrain = createStrainFromFile(node)
            if newStrain != None:
                print('Successfully created the following strain from data file: %s' % (newStrain.getName()))
                strains.append(newStrain)
                return newStrain

    if leftSibling != None and leftSibling.getSequence() != None and len(leftSibling.getSequence()) > 0 and rightSibling != None and rightSibling.getSequence() != None and len(rightSibling.getSequence()) > 0:
        print('Strains compared: %s, %s' % (leftSibling.getName(), rightSibling.getName()))

        neighborStrain = None
        if parentNode != None:
            node.name = 'Processing'
            #Helps determine which side we need to traverse
            if len(parentNode.clades) > 0 and parentNode.clades[0].name != 'Processing':
                neighborStrain = getNeighborStrain(parentNode.clades[0])
            elif len(parentNode.clades) > 1 and parentNode.clades[1].name != 'Processing':
                neighborStrain = getNeighborStrain(parentNode.clades[1])
            #Put the name back the way it was
            node.name = None

        if neighborStrain != None:
            print('The following neighbor will be used during the comparison: %s' % (neighborStrain.getName()))
        else:
            print('No neighbor found!')
            
        ancestralOperons, events = processStrains(leftSibling, rightSibling, neighborStrain)
        globals.ancestralCounter += 1
        node.name = 'Ancestor %s' % (globals.ancestralCounter)
        
        ancestor = Strain('Ancestor %s' % (globals.ancestralCounter), ancestralOperons, [leftSibling.getName(), rightSibling.getName()], [], {})
        ancestor.setEvents(events)
        strains.append(ancestor)

        return ancestor

    #If the left child has a sequence, return it
    elif not(leftSibling == None) and len(leftSibling.getSequence()) > 0:
        return leftSibling

    #If the right child has a sequence, return it
    elif not(rightSibling == None) and len(rightSibling.getSequence()) > 0:
        return rightSibling

    #If neither has a sequence, return None
    else:
        return None

######################################################
# getNeighborStrain
# Parameters:
# Description: Finds the first available strain with data
######################################################
def getNeighborStrain(currNode):
    global strains
    neighbor = None

    if neighbor is None and currNode.name != None and len(currNode.name) > 0:
        #Retrieves strain based on the name if it exists in the list
        filteredList = iter(filter(lambda x: x.name == currNode.name, strains))
        neighbor = next(filteredList, None)

        #Create strain from data file if not an Ancestor (avoids reading in and processing files a second time in the future)
        if (neighbor == None and not('Ancestor' in currNode.name)):
            neighbor = createStrainFromFile(currNode)
            if neighbor != None:
                print('Successfully created new strain from data file while looking for neighbor: %s' % (neighbor.getName()))
                strains.append(neighbor)

    if neighbor == None and len(currNode.clades) > 0:
        neighbor = getNeighborStrain(currNode.clades[0])

    if neighbor == None and len(currNode.clades) > 1:
        neighbor = getNeighborStrain(currNode.clades[1])

    return neighbor

######################################################
# createStrainFromFile
# Parameters: node - node from Newick tree to create the strain from
# Description: Reads-in data file for strain and returns a Strain object
######################################################
def createStrainFromFile(node):
    strain = None

    if os.path.isdir(node.name):
        if os.path.isfile(node.name + '/sequence.txt'):
            fileGeneSequence = open(node.name + '/sequence.txt', 'r').read()
            operons, operonPositions, singletonDict, allGenes = processFileSequence(fileGeneSequence)

            #Random newline characters in file are causing the alignment scores to be messed up
            for x in range(0, len(operons)):
                operons[x] = (operons[x]).replace('\r', '').replace('\n', '')
            for x in range(0, len(allGenes)):
                allGenes[x] = (allGenes[x]).replace('\r', '').replace('\n', '')
            for key in singletonDict.keys():
                singletonDict[key] = (singletonDict[key]).replace('\r', '').replace('\n', '')

            strain = Strain(node.name, operons, [], operonPositions, singletonDict)
            strain.setGenes(allGenes)
        else:
            print('No sequence file found for node: %s' % node.name)
    else:
        print('No directory found for node: %s' % node.name)

    return strain

######################################################
# processFileSequence
# Parameters: sequence - string sequence of strain
# Description: Processes a string sequence passed in into any array
######################################################
def processFileSequence(sequence):
    geneList = []
    operonList = []
    singletonList = {}
    index = 0
    geneIndex = 0
    operonPositions = []

    while index < len(sequence):

        #get the origin and termus (remove if not needed)
        if sequence[index] == '<':
            startIndex = index

            while sequence[index] != '>':
                index += 1

            #increment the index to include the >
            index += 1
            # operonList.append(sequence[startIndex:index])
            #decrement the geneIndex to not include <> values
            geneIndex -= 1

        #Operon
        elif (sequence[index] == '[') or (sequence[index] == '-' and sequence[index + 1] == '['):
            startIndex = index
            operonPositions.append(geneIndex)

            while sequence[index] != ']':
                if sequence[index] == ',':
                    geneIndex += 1
                index += 1

            #increment the index to include the ]
            index += 1
            operonList.append(sequence[startIndex:index])
            if sequence[startIndex] == '[':
                geneList.extend([gene.strip() for gene in sequence[startIndex+1:index-1].split(',')])
            else:
                geneList.extend([gene.strip() for gene in sequence[startIndex+2:index-1].split(',')])

        #Singleton
        #TODO need to verify this code (what happens if singleton first instead of an operon?)
        if index < len(sequence) and sequence[index] == ',':
            geneIndex += 1
            if sequence[index+2] != '[' and sequence[index+2] != '<' and sequence[index+3] != '[':
                index += 1
                startIndex = index

                while sequence[index] != ',' and index < len(sequence)-1:
                    index += 1

                #Add singleton gene
                operonList.append(sequence[startIndex:index])
                operonPositions.append(geneIndex)

                singletonList[str(geneIndex)] = sequence[startIndex:index]
                geneList.append(sequence[startIndex:index].strip().replace("-", ""))
                index -= 1
        index += 1

    return operonList, operonPositions, singletonList, geneList

######################################################
#                       main
######################################################
globals.initialize()
createFile(outputFile)

startTime = time.time()

print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree recursively reconstructing ancestral genomes
result = traverseNewickTree(newickTree.clade, None)

endTime = time.time()
print('Total number of seconds: %s\nProcessing Complete!' % (endTime - startTime))