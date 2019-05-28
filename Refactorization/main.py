import time
import copy
import os.path
import globals
from Bio import Phylo
from FileService import createFile
from FileService import processSequence
from LineageSummary import LineageSummary
from SequenceService import createDotPlot
from FileService import outputTotalsToFile
from SequenceService import createBarGraph
from BacterialStrain import BacterialStrain
from FragmentService import determineRegions
from FileService import outputStrainDetailsToFile
from SequenceService import normalizeIndexesForDotPlot
from SequenceService import updateGlobalDeletionCounter
from SequenceService import updateGlobalDuplicationCounter
from LocalAlignmentModule import findOrthologsByLocalAlignment
from GlobalAlignmentModule import findOrthologsByGlobalAlignment
from SelfGlobalAlignmentModule import findOrthologsBySelfGlobalAlignment
from SequenceService import updateGlobalInversionSizeDistributionCounter
from SequenceService import updateGlobalTranspositionSizeDistributionCounter
from FragmentService import determineAncestralFragmentArrangementUsingNeighbor
from FragmentService import determineAncestralFragmentArrangementWithoutNeighbor
from SequenceService import updateGlobalInvertedTranspositionSizeDistributionCounter

#Application parameters
newickFileName = 'Bacillus_Tree.dnd' #Name of newick tree file

#Global vaiables used by the script
strains = [] #Global variable that stores all of the strains in the phylogenetic tree
outputFileName = 'ApplicationOutput.txt' #Name of output file

################################################################
################# Main Fuctions ################################
################################################################

######################################################
# createAncestor
# Parameters:
# Description: Takes two related strains and a close neighbor and constructs the ancestral node
######################################################
def createAncestor(strain1, strain2, neighborStrain):
    globals.ancestralCounter+=1
    ancestor = None
    ancestralName = 'Ancestor ' + str(globals.ancestralCounter)
    ancestralFragments = None

    strain1Copy = copy.deepcopy(strain1) #Do a deep copy of object for when we compare to the neighbor
    neighborCopy = copy.deepcopy(neighborStrain) #Do a deep copy of the neighbor as well b/c we don't want to store those comparisons in the strain either

    print('Performing a series of alignments for the following strains: %s, %s' % (strain1.name, strain2.name))
    events, duplicatesStrain1, duplicatesStrain2 = constructEvents(strain1, strain2)

    print('Constructing dot plot for the following strains: %s, %s' % (strain1.name, strain2.name))
    points, lostPoints = normalizeIndexesForDotPlot(events, duplicatesStrain1, duplicatesStrain2, strain1, strain2)
    createDotPlot(points, strain1, strain2)

    createBarGraph(strain1.duplicationCounts, 'Distribution of Duplications for %s'%(strain1.name))
    createBarGraph(strain2.duplicationCounts, 'Distribution of Duplications for %s'%(strain2.name))
    createBarGraph(strain1.deletionCounts, 'Distribution of Deletions for %s'%(strain1.name)) #Remember! Deletions refer to the other strain!
    createBarGraph(strain2.deletionCounts, 'Distribution of Deletions for %s'%(strain2.name)) #Remember! Deletions refer to the other strain!

    #Compute and output the inverted, transposed, and inverted transposed regions
    FCR, TR, IR, ITR = determineRegions(points)
    #FCR, TR, IR, ITR, LR = computeOperonArrangements(events)  OLD VERSION

    #inversionDetails1, inversionDetails2 = computeRegionDetails(IR, 'Inversion:')
    #transpositionDetails1, transpositionDetails2 = computeRegionDetails(TR, 'Transposition:')
    #invertedTransposedDetails1, invertedTransposedDetails2 = computeRegionDetails(ITR, 'Inverted Transposition:')

    #TODO singletons should not have brackets? Still deciding but this should not affect anything

    #Compare one of the siblings to the neighbor if one exists
    if neighborCopy != None:
        print('Now performing a series of alignments between the nighboring strains: %s, %s' % (strain1Copy.name, neighborCopy.name))
        neighborEvents, duplicatesStrain1Copy, duplicatesStrainNeighbor = constructEvents(strain1Copy, neighborCopy)

        print('Constructing dot plot for the neighboring strains: %s, %s' % (strain1Copy.name, neighborCopy.name))
        neighborPoints, neighborLostPoints = normalizeIndexesForDotPlot(neighborEvents, duplicatesStrain1Copy, duplicatesStrainNeighbor, strain1Copy, neighborCopy)
        #createDotPlot(neighborPoints, strain1Copy, neighborCopy)

        #Compute the various regions for the neighbor
        #NFCR, NTR, NIR, NITR, NLR = computeOperonArrangements(neighborEvents) OLD VERSION
        NFCR, NTR, NIR, NITR = determineRegions(neighborPoints)
        ancestralFragments, strain1, strain2 = determineAncestralFragmentArrangementUsingNeighbor(FCR, TR, IR, ITR, lostPoints, NFCR, NTR, NIR, NITR, neighborLostPoints, strain1, strain2)
    else:
        if neighborCopy == None:
            print('No neighbor found!')
        elif len(TR) == 0 and len(IR) == 0 or len(ITR) == 0:
            print('No inverted or transposed regions detected!!')
        ancestralFragments, strain2 = determineAncestralFragmentArrangementWithoutNeighbor(FCR, TR, IR, ITR, lostPoints, strain2)
    
    #Computes the total number of inversions, transpositions, inverted transpositions
    globals.inversionCounter += len(IR)
    globals.transposedCounter += len(TR)
    globals.invertedTransposedCounter += len(ITR)
    
    #Increments the counters for the size distributions for each event type
    updateGlobalDeletionCounter(strain1)
    updateGlobalDeletionCounter(strain2)
    updateGlobalDuplicationCounter(strain1)
    updateGlobalDuplicationCounter(strain2)
    updateGlobalInversionSizeDistributionCounter(strain1)
    updateGlobalInversionSizeDistributionCounter(strain2)
    updateGlobalTranspositionSizeDistributionCounter(strain1)
    updateGlobalTranspositionSizeDistributionCounter(strain2)
    updateGlobalInvertedTranspositionSizeDistributionCounter(strain1)
    updateGlobalInvertedTranspositionSizeDistributionCounter(strain2)
    
    #Append all details to file here
    outputStrainDetailsToFile(outputFileName, strain1)
    outputStrainDetailsToFile(outputFileName, strain2)
    
    ancestor = BacterialStrain(ancestralName, ancestralFragments)
    return ancestor

######################################################
# initializeTracker
# Parameters: Strain
# Description: Initializes an array to False for all genome fragments which will be used to track which fragments have been marked
######################################################
def initializeTracker(strain):
    coverageTracker = {}
    for x in range(0, len(strain.genomeFragments)):
        coverageTracker[x] = False
    return coverageTracker

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
# Description: Constructs the orthologous events between two provided strains
######################################################
def constructEvents(strain1, strain2):

    events = [] #Stores all of the orthologous events also keeps track of deletions and duplications
    duplicationEvents1 = [] #Stores the duplication events for strain 1
    duplicationEvents2 = [] #Stores the duplication events for strain 2
    coverageTracker1 = initializeTracker(strain1) #Tracks which operons have been marked off in strain 1
    coverageTracker2 = initializeTracker(strain2) #Tracks which operons have been marked off in strain 2

    print('Performing global alignment with: %s, %s' % (strain1.name, strain2.name))
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)

    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))

    #Local Alignment operation
    if numRemainingOperons1 > 0 and numRemainingOperons2 > 0:
        print('Performing local alignment with: %s, %s' % (strain1.name, strain2.name))
        localAlignmentEvents, coverageTracker1, coverageTracker2, localAlignmentCounter, strain1, strain2 = findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
        print('Number of orthologous operons identified using Local Alignment %s' % (localAlignmentCounter))

        numRemainingOperons1 = countRemainingOperons(coverageTracker1)
        numRemainingOperons2 = countRemainingOperons(coverageTracker2)
        print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))
        if len(localAlignmentEvents) > 0:
            events.extend(localAlignmentEvents)

    #Self Global Alignment
    if numRemainingOperons1 > 0:
        #Remember to insert the deletions into the sibling (that's how we defined it)
        duplicationEvents1, lossEvents1, coverageTracker1, strain1, strain2 = findOrthologsBySelfGlobalAlignment(strain1, coverageTracker1, strain2)
        print('%s, duplicates identified %s and losses identified %s' % (strain1.name, len(duplicationEvents1), len(lossEvents1)))
        if len(lossEvents1) > 0:
            events.extend(lossEvents1)

    if numRemainingOperons2 > 0:
        #Remember to insert the deletions into the sibling (that's how we defined it)
        duplicationEvents2, lossEvents2, coverageTracker2, strain2, strain1 = findOrthologsBySelfGlobalAlignment(strain2, coverageTracker2, strain1)
        print('%s, duplicates identified %s and losses identified %s' % (strain2.name, len(duplicationEvents2), len(lossEvents2)))
        if len(lossEvents2) > 0:
            events.extend(lossEvents2)

    #Verify there's no unmarked operons at this point
    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    if numRemainingOperons1 > 0 or numRemainingOperons2 > 0:
        print('Error! There are unmarked operons remaining!')

    return events, duplicationEvents1, duplicationEvents2
######################################################
# getNeighborStrain
# Parameters:
# Description: Finds the first available strain with data
######################################################
def getNeighborStrain(currNode):
    global strains
    neighbor = None

    if neighbor == None and currNode.name != None and len(currNode.name) > 0:
        #Retrieves strain based on the name if it exists in the list
        filteredList = iter(filter(lambda x: x.name == currNode.name, strains))
        neighbor = next(filteredList, None)

        #Create strain from data file if not an Ancestor (avoids reading in and processing files a second time in the future)
        if (neighbor == None and not('Ancestor' in currNode.name)):
            neighbor = createStrainFromFile(currNode)
            if neighbor != None:
                print('Successfully created new strain from data file while looking for neighbor: %s' % (neighbor.name))
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
            genome = open(node.name + '/sequence.txt', 'r').read()
            strain = processSequence(node.name, genome)
        else:
            print('No sequence file found for node: %s' % node.name)
    else:
        print('No directory found for node: %s' % node.name)

    return strain

######################################################
# computeLineageCost
# Parameters:
# Description: Traverses newick tree in post order to compute the cost of a lineage
######################################################
def computeLineageCost(node, targetName, lineageCost):
    newLineageCost = LineageSummary(targetName)
    
    if lineageCost != None: #Insert the previous costs
        newLineageCost.totalCodonMismatches = lineageCost.totalCodonMismatches
        newLineageCost.totalSubstitutions = lineageCost.totalSubstitutions
        
    #Now add in the costs for the current node
    filteredList = iter(filter(lambda x: x.name == node.name, strains))
    foundStrain = next(filteredList, None)
    if foundStrain != None:
        count = foundStrain.codonMismatchDetails.count(';')
        newLineageCost.totalCodonMismatches += count
        count = foundStrain.substitutionDetails.count(';')
        newLineageCost.totalSubstitutions += count
    else:
        print('Error! Unable to find the following strain: %s' % (node.name))
        return None
    
    if node.name == targetName: #Check if this is our target
        print('Found  node in newick tree! The found node is: %s' % (targetName))
        return newLineageCost
    
    if len(node.clades) > 0:
        temp = computeLineageCost(node.clades[0], targetName, newLineageCost)
        if temp != None:
            return temp #If we found our target
        else:
            if len(node.clades) > 1:
                temp = computeLineageCost(node.clades[1], targetName, newLineageCost)
                if temp != None:
                    return temp #If we found our target
    return None
    
######################################################
# traverseNewickTree
# Parameters: node - Strain being currently processed, parentNode - direct ancestor of node
# Description: Traverses a provided newick tree in post order traversal
######################################################
def traverseNewickTree(node, parentNode):
    global strains

    leftSibling = None
    rightSibling = None

    #Check there's any descendants
    if len(node.clades) > 0:
        leftSibling = traverseNewickTree(node.clades[0], node)
        if len(node.clades) > 1:
            rightSibling = traverseNewickTree(node.clades[1], node)

    #Retrieve strain or create one from the data file
    if node.name != None and len(node.name) > 0:
        #This code attempts to retrieve the strain from the strain list
        filteredList = iter(filter(lambda x: x.name == node.name, strains))
        foundStrain = next(filteredList, None)

        if (foundStrain != None):
            print('Retrieving strain from strains list: %s' % (foundStrain.name))
            return foundStrain
        else:
            newStrain = createStrainFromFile(node)
            if newStrain != None:
                print('Successfully created the following strain from data file: %s' % (newStrain.name))
                strains.append(newStrain)
                return newStrain

    #Case 1: Both siblings exist therefore we need to construct their ancestor
    if leftSibling != None and rightSibling != None and leftSibling.genomeFragments != None and len(leftSibling.genomeFragments) > 0 and rightSibling.genomeFragments != None and len(rightSibling.genomeFragments) > 0:
        print('The following siblings will be compared, %s, %s...' % (leftSibling.name, rightSibling.name))

        neighborStrain = None #Neighbor strain
        if parentNode != None:
            node.name = 'Processing' #Helps determine whether we go left or right to get the neighbor
            if len(parentNode.clades) > 0 and parentNode.clades[0].name != 'Processing':
                neighborStrain = getNeighborStrain(parentNode.clades[0])
            elif len(parentNode.clades) > 1 and parentNode.clades[1].name != 'Processing':
                neighborStrain = getNeighborStrain(parentNode.clades[1])
            node.name = None #Put the name back the way it was so we don't mess up anything

        if neighborStrain != None:
            print('In addition the following neighbor will be used during the comparison, %s' % (neighborStrain.name))
        else:
            print('No neighbor found!')

        ancestor = createAncestor(leftSibling, rightSibling, neighborStrain)
        strains.append(ancestor)

        return ancestor

    #Case 2: Only the left sibling exists so return it
    elif leftSibling != None and leftSibling.genomeFragments != None and len(leftSibling.genomeFragments) > 0:
        return leftSibling
    #Case 3: Only the right sibling exists so return it
    elif rightSibling != None and rightSibling.genomeFragments != None and len(rightSibling.genomeFragments) > 0:
        return rightSibling
    #Case 4: None of the siblings exist so return NULL
    else:
        return None

######################################################
#                       main
######################################################
print('Starting application...')
startTime = time.time()

print('Reading newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

globals.initialize() #Initialize the globals file
createFile(outputFileName, newickTree) #Creates file where data will be output

#Traverses the newick tree recursively reconstructing ancestral genomes
print('Traversing newick tree...')
result = traverseNewickTree(newickTree.clade, None)

#Output the totals for the computation to console and file
outputTotalsToFile(outputFileName)

#Output newick tree after the ancestors have been added to it
Phylo.draw(newickTree)

endTime = time.time()
totalTime = endTime - startTime
print('Total time (in seconds): %s' % (totalTime))

#TODO compute lineage
target = 'NC_014019'
print('Computing lineage cost for: %s' % (target))
lineageCost = computeLineageCost(newickTree.clade, target, None)
if lineageCost != None:
    print('Successfully found and computed the lineage for: %s' % (target))

print('Ending application...')