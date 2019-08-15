import time
import copy
import os.path
import sys
import globals
from Bio import Phylo
from FileService import createFile
from FileService import processSequence
from LineageSummary import LineageSummary
from SequenceService import createDotPlot
from FileService import outputTotalsToFile
from SequenceService import createBarGraph
from FileService import outputGenomeToFile
from BacterialStrain import BacterialStrain
from FragmentService import determineRegions
from FileService import outputStrainDetailsToFile
from SequenceService import normalizeIndexesForDotPlot
from SequenceService import updateGlobalDeletionCounter
from GlobalAlignmentModule import reduceSingletonDeletions
from SequenceService import updateGlobalDuplicationCounter
from SequenceService import updateGlobalSubstitutionCounter
from SequenceService import updateGlobalCodonMismatchCounter
from LocalAlignmentModule import findOrthologsByLocalAlignment
from GlobalAlignmentModule import findOrthologsByGlobalAlignment
from SelfGlobalAlignmentModule import findOrthologsBySelfGlobalAlignment
from SequenceService import updateGlobalInversionSizeDistributionCounter
from SequenceService import updateGlobalTranspositionSizeDistributionCounter
from FragmentService import determineAncestralFragmentArrangementUsingNeighbor
from FragmentService import determineAncestralFragmentArrangementWithoutNeighbor
from SequenceService import updateGlobalInvertedTranspositionSizeDistributionCounter

#Application parameters
newickFileName = 'tree.dnd' #Name of newick tree file

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

    if globals.printToConsole:
        print('Performing a series of alignments for the following strains: %s, %s' % (strain1.name, strain2.name))
        
    globals.enableDeletionReversions = True #Only do the backtrace between these two strains!
    globals.enableSelfAlignmentDetails = True
    
    events, duplicatesStrain1, duplicatesStrain2 = constructEvents(strain1, strain2)
    
    globals.enableSelfAlignmentDetails = False
    globals.enableDeletionReversions = False

    if globals.printToConsole:
        print('Constructing dot plot for the following strains: %s, %s' % (strain1.name, strain2.name))
    points, lostPoints = normalizeIndexesForDotPlot(events, duplicatesStrain1, duplicatesStrain2, strain1, strain2)
    if globals.printToConsole:
        createDotPlot(points, strain1, strain2, testFileName)
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

    #Compare one of the siblings to the neighbor if one exists
    if neighborCopy != None:
        if globals.printToConsole:
            print('Now performing a series of alignments between the nighboring strains: %s, %s' % (strain1Copy.name, neighborCopy.name))
        neighborEvents, duplicatesStrain1Copy, duplicatesStrainNeighbor = constructEvents(strain1Copy, neighborCopy)
        if globals.printToConsole:
            print('Constructing dot plot for the neighboring strains: %s, %s' % (strain1Copy.name, neighborCopy.name))
        neighborPoints, neighborLostPoints = normalizeIndexesForDotPlot(neighborEvents, duplicatesStrain1Copy, duplicatesStrainNeighbor, strain1Copy, neighborCopy)
        #createDotPlot(neighborPoints, strain1Copy, neighborCopy)

        #Compute the various regions for the neighbor
        #NFCR, NTR, NIR, NITR, NLR = computeOperonArrangements(neighborEvents) OLD VERSION
        NFCR, NTR, NIR, NITR = determineRegions(neighborPoints)
        ancestralFragments, strain1, strain2 = determineAncestralFragmentArrangementUsingNeighbor(FCR, TR, IR, ITR, lostPoints, NFCR, NTR, NIR, NITR, neighborLostPoints, strain1, strain2)
    else:
        if neighborCopy == None:
            if globals.printToConsole:
                print('No neighbor found!')
        elif len(TR) == 0 and len(IR) == 0 or len(ITR) == 0:
            if globals.printToConsole:
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

    #Increment counters (only need to do the count only once otherwise it leads to double counts ie x2 number of events)
    #updateGlobalCodonMismatchCounter(strain1)
    updateGlobalCodonMismatchCounter(strain2)
    #updateGlobalSubstitutionCounter(strain1)
    updateGlobalSubstitutionCounter(strain2)

    #Append all details to file here
    #outputStrainDetailsToFile(outputFileName, strain1)
    #outputStrainDetailsToFile(outputFileName, strain2)

    ancestor = BacterialStrain(ancestralName, ancestralFragments)

    if globals.printToConsole:
        print(strain1.name)
        for frag in strain1.genomeFragments:
            print(frag.originalSequence)
        print(strain2.name)
        for frag in strain2.genomeFragments:
            print(frag.originalSequence)

    ####################################
    #Handle the Codon Mismatches here##
    ###################################
    if '#' in strain1.codonMismatchDetails:
        newDetails1 = 'Codon Mismatch:'
        newDetails2 = 'Codon Mismatch:'

        line1 = strain1.codonMismatchDetails.replace('Codon Mismatch:', '').strip()
        line2 = strain2.codonMismatchDetails.replace('Codon Mismatch:', '').strip()

        subsList1 = filter(None, line1.split(';')) #Ensures we don't have a list with an empty string as an element
        subsList2 = filter(None, line2.split(';'))

        #For each substitution in the list
        for w in range(0, len(subsList1)):
            gene1, idNumber1, position1 = parseDetails(subsList1[w])
            gene2, idNumber2, position2 = parseDetails(subsList2[w])
            processed = False #Tracks whether the current codon mismatch was handled

            #Check if we have a neighbor
            if neighborCopy:
                #Check if the same codon mismatch occurred when comparing to the neighbor
                if '#' in strain1Copy.codonMismatchDetails:
                    line3 = strain1Copy.codonMismatchDetails.replace('Codon Mismatch:', '').strip()
                    subsList3 = filter(None, line3.split(';'))
                    for v in range(0, len(subsList3)):
                        gene3, idNumber3, position3 = parseDetails(subsList3[v])
                        if gene1 == gene3 and position1 == position3:
                            #We found the same codon mismatch when comparing with the neighbor, therefore we should keep strain 2's verison of the gene!
                            processed = True
                            fragments = ancestor.genomeFragments
                            for fragment in fragments:
                                if idNumber1 in fragment.originalSequence:
                                    fragment.originalSequence = fragment.originalSequence.replace(gene1 + '-' + idNumber1, gene2) #Put in strain 2's gene
                                    for m in range(0, len(fragment.sequence)):
                                        if idNumber1 in fragment.sequence[m]:
                                            fragment.sequence[m] = gene2
                                            break
                                    break
            if processed:
                #We found the codon mismatch and swapped with strain 2's gene therefore strain 1's gene was the codon mismatch so put the codon mismatch details in strain1
                newDetails1+= gene1 + ' ' + position1 + ';'
            else:
                #We were not able to find the same codon mismatch either due to there being no neighbor or it was just not there. So just assume strain 2 is the codon mismatch
                newDetails2+= gene2 + ' ' + position2 + ';'
                fragments = ancestor.genomeFragments
                for fragment in fragments:
                    if idNumber1 in fragment.originalSequence:
                        fragment.originalSequence = fragment.originalSequence.replace(gene1 + '-' + idNumber1, gene1) #Put in strain 1's gene
                        for m in range(0, len(fragment.sequence)):
                            if idNumber1 in fragment.sequence[m]:
                                fragment.sequence[m] = gene1
                                break
                        break
        #Insert the new details about the substitution
        strain1.codonMismatchDetails = newDetails1
        strain2.codonMismatchDetails = newDetails2

    ################################
    #Handle the substitutions here##
    ################################
    if '@' in strain1.substitutionDetails:
        newDetails1 = 'Substitution:'
        newDetails2 = 'Substitution:'

        line1 = strain1.substitutionDetails.replace('Substitution:', '').strip()
        line2 = strain2.substitutionDetails.replace('Substitution:', '').strip()

        subsList1 = filter(None, line1.split(';')) #Ensures we don't have a list with an empty string as an element
        subsList2 = filter(None, line2.split(';'))

        #For each substitution in the list
        for w in range(0, len(subsList1)):
            gene1, idNumber1, position1 = parseDetails(subsList1[w])
            gene2, idNumber2, position2 = parseDetails(subsList2[w])
            processed = False #Tracks whether the current substitution was handled

            #Check if we have a neighbor
            if neighborCopy:
                #Check if the same substitution occurred when comparing to the neighbor
                if '@' in strain1Copy.substitutionDetails:
                    line3 = strain1Copy.substitutionDetails.replace('Substitution:', '').strip()
                    subsList3 = filter(None, line3.split(';'))
                    for v in range(0, len(subsList3)):
                        gene3, idNumber3, position3 = parseDetails(subsList3[v])
                        if gene1 == gene3 and position1 == position3:
                            #We found the same substitution when comparing with the neighbor, therefore we should keep strain 2's verison of the gene!
                            processed = True
                            fragments = ancestor.genomeFragments
                            for fragment in fragments:
                                if idNumber1 in fragment.originalSequence:
                                    fragment.originalSequence = fragment.originalSequence.replace(gene1 + '-' + idNumber1, gene2) #Put in strain 2's gene
                                    for m in range(0, len(fragment.sequence)):
                                        if idNumber1 in fragment.sequence[m]:
                                            fragment.sequence[m] = gene2
                                            break
                                    break
            if processed:
                #We found the substitution and swapped with strain 2's gene therefore strain 1's gene was the substituion so put the substitution details in strain1
                newDetails1+= gene1 + ' ' + position1 + ';'
            else:
                #We were not able to find the same substitution either due to there being no neighbor or it was just not there. So just assume strain 2 is the substitution
                newDetails2+= gene2 + ' ' + position2 + ';'
                fragments = ancestor.genomeFragments
                for fragment in fragments:
                    if idNumber1 in fragment.originalSequence:
                        fragment.originalSequence = fragment.originalSequence.replace(gene1 + '-' + idNumber1, gene1) #Put in strain 1's gene
                        for m in range(0, len(fragment.sequence)):
                            if idNumber1 in fragment.sequence[m]:
                                fragment.sequence[m] = gene1
                                break
                        break
        #Insert the new details about the substitution
        strain1.substitutionDetails = newDetails1
        strain2.substitutionDetails = newDetails2

    
    #Add any codon mismatches from the self global alignment as those details were stored in another variable so it doesn't mess with codon mismatches and substitution handlers in the previous 2 for loops
    strain1.codonMismatchDetails += strain1.tempCodonDetails
    strain2.codonMismatchDetails += strain2.tempCodonDetails
    strain1.substitutionDetails += strain1.tempSubstitutionDetails
    strain2.substitutionDetails += strain2.tempSubstitutionDetails
    
    return ancestor

######################################################
# parseDetails
# Parameters:
# Description: Parses the gene, unique Id number, and position from a given substitution
######################################################
def parseDetails(line):
    temp = line.split('-')
    gene = temp[0].strip() #Extracts the substituted gene

    temp = temp[1].split(' ')
    idNumber = temp[0] #Extracts unique Id number
    position = temp[1] #Extracts the position

    return gene, idNumber, position

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
    lossEvents1 = []
    lossEvents2 = []
    newEvents = []

    if globals.printToConsole:
        print('Performing global alignment with: %s, %s' % (strain1.name, strain2.name))
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)

    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    if globals.printToConsole:
        print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))

    #Local Alignment operation
#    if numRemainingOperons1 > 0 and numRemainingOperons2 > 0:
#        if globals.printToConsole:
#            print('Performing local alignment with: %s, %s' % (strain1.name, strain2.name))
#        localAlignmentEvents, coverageTracker1, coverageTracker2, localAlignmentCounter, strain1, strain2 = findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
#        if globals.printToConsole:
#            print('Number of orthologous operons identified using Local Alignment %s' % (localAlignmentCounter))
#
#        numRemainingOperons1 = countRemainingOperons(coverageTracker1)
#        numRemainingOperons2 = countRemainingOperons(coverageTracker2)
#        if globals.printToConsole:
#            print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))
#        if len(localAlignmentEvents) > 0:
#            events.extend(localAlignmentEvents)

    #Self Global Alignment
    if numRemainingOperons1 > 0:
        #Remember to insert the deletions into the sibling (that's how we defined it)
        duplicationEvents1, lossEvents1, coverageTracker1, strain1, strain2 = findOrthologsBySelfGlobalAlignment(strain1, coverageTracker1, strain2)
        if globals.printToConsole:
            print('%s, duplicates identified %s and losses identified %s' % (strain1.name, len(duplicationEvents1), len(lossEvents1)))

    if numRemainingOperons2 > 0:
        #Remember to insert the deletions into the sibling (that's how we defined it)
        duplicationEvents2, lossEvents2, coverageTracker2, strain2, strain1 = findOrthologsBySelfGlobalAlignment(strain2, coverageTracker2, strain1)
        if globals.printToConsole:
            print('%s, duplicates identified %s and losses identified %s' % (strain2.name, len(duplicationEvents2), len(lossEvents2)))

    #Try reducing the number of singleton deletions
#    if len(lossEvents1) > 0 or len(lossEvents2) > 0:
#        lossEvents1, lossEvents2, newEvents = reduceSingletonDeletions(lossEvents1, lossEvents2, coverageTracker1, coverageTracker2, strain1, strain2)
#    if len(newEvents) > 0:
#        events.extend(newEvents)
    if len(lossEvents1) > 0:
        events.extend(lossEvents1)
    if len(lossEvents2) > 0:
        events.extend(lossEvents2)

    #Verify there's no unmarked operons at this point
    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    if numRemainingOperons1 > 0 or numRemainingOperons2 > 0:
        if globals.printToConsole:
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
                if globals.printToConsole:
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

    if os.path.isdir(testFileName + node.name):
        if os.path.isfile(testFileName + node.name + '/sequence.txt'):
            genome = open(testFileName + node.name + '/sequence.txt', 'r').read()
            strain = processSequence(node.name, genome)
        else:
            if globals.printToConsole:
                print('No sequence file found for node: %s' % node.name)
    else:
        if globals.printToConsole:
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

    if node.name != None: #Sometimes there's a none type. not sure why
        #Now add in the costs for the current node
        filteredList = iter(filter(lambda x: x.name == node.name, strains))
        foundStrain = next(filteredList, None)

        if foundStrain != None:
            count = foundStrain.codonMismatchDetails.count(';')
            newLineageCost.totalCodonMismatches += count
            count = foundStrain.substitutionDetails.count(';')
            newLineageCost.totalSubstitutions += count
        if node.name == targetName: #Check if this is our target
            if globals.printToConsole:
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
# traverseNewickTreeAndOutputToFile
# Parameters: node - Strain being currently processed
# Description: Traverses a provided newick tree in post order traversal to output the appropriate details to the output file
######################################################
def traverseNewickTreeAndOutputToFile(node):
    if len(node.clades) > 0:
        traverseNewickTreeAndOutputToFile(node.clades[0])
        if len(node.clades) > 1:
            traverseNewickTreeAndOutputToFile(node.clades[1])
    if node.name != None and len(node.name) > 0:
        filteredList = iter(filter(lambda x: x.name == node.name, strains))
        foundStrain = next(filteredList, None)
        if foundStrain != None:
            outputStrainDetailsToFile(outputFileName, foundStrain)
            outputGenomeToFile(node.name + ".txt", foundStrain)

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
            if globals.printToConsole:
                print('Retrieving strain from strains list: %s' % (foundStrain.name))
            return foundStrain
        else:
            newStrain = createStrainFromFile(node)
            if newStrain != None:
                if globals.printToConsole:
                    print('Successfully created the following strain from data file: %s' % (newStrain.name))
                strains.append(newStrain)
                return newStrain

    #Case 1: Both siblings exist therefore we need to construct their ancestor
    if leftSibling != None and rightSibling != None and leftSibling.genomeFragments != None and len(leftSibling.genomeFragments) > 0 and rightSibling.genomeFragments != None and len(rightSibling.genomeFragments) > 0:
        if globals.printToConsole:
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
            if globals.printToConsole:
                print('In addition the following neighbor will be used during the comparison, %s' % (neighborStrain.name))
        else:
            if globals.printToConsole:
                print('No neighbor found!')

        ancestor = createAncestor(leftSibling, rightSibling, neighborStrain)
        node.name = ancestor.name
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
def main():
    globals.initialize() #Initialize the globals file
    
    global newickFileName
    global outputFileName
    global testFileName
    
    if len(sys.argv) != 3:
        print "WARNING: Must provide a Newick tree and test folder name. Exiting..."
        sys.exit(0)
    
    newickFileName = sys.argv[1]
    if newickFileName == "tree2LeafNeighbour.dnd":
        outputFileName = sys.argv[2] + "/ApplicationNeighbourOutput.txt"
    else:
        outputFileName = sys.argv[2] + "/ApplicationOutput.txt"
    testFileName = sys.argv[2] + '/'
    
    print('Starting application...')
    startTime = time.time()
    
    if globals.printToConsole:
        print('Reading newick tree from file: %s...' % (newickFileName))
    newickTree = Phylo.read(newickFileName, 'newick')
    if globals.printToConsole:
        Phylo.draw(newickTree)
    
    
    globals.strains = strains #Assign pointer to the global strains array so we can access it anywhere
    createFile(outputFileName, newickTree) #Creates file where data will be output
    
    #Traverses the newick tree recursively reconstructing ancestral genomes
    if globals.printToConsole:
        print('Traversing newick tree...')
    result = traverseNewickTree(newickTree.clade, None)
    
    endTime = time.time()
    totalTime = endTime - startTime
    
    #Output ancestral genome to console
    if globals.printToConsole:
        print('This is the root ancestral genome!')
        
    root = newickTree.clade
    rootGenome = []
    if newickFileName == "tree2LeafNeighbour.dnd":
        if len(root.clades) == 2:
            child = root.clades[0]
            if len(child.clades) != 2:
                child = root.clades[1]
                neighbour = root.clades[0]
            else:
                neighbour = root.clades[1]
            if child.name != None and len(child.name) > 0:
                filteredList = iter(filter(lambda x: x.name == child.name, strains))
                foundStrain = next(filteredList, None)
                if foundStrain != None:
                    ancestralFragments = foundStrain.genomeFragments
                    rootGenome = ', '.join(fragment.originalSequence for fragment in ancestralFragments)
                        
            with open(testFileName + "appNeighbourRoot.txt", "w+") as f:
                f.write(rootGenome)
            neighbour.name = ''
            child.name = ''
    else:
        if root.name != None and len(root.name) > 0:
            filteredList = iter(filter(lambda x: x.name == root.name, strains))
            foundStrain = next(filteredList, None)
            if foundStrain != None:
                ancestralFragments = foundStrain.genomeFragments
                rootGenome = ', '.join(fragment.originalSequence for fragment in ancestralFragments)
                    
        with open(testFileName + "appRoot.txt", "w+") as f:
            f.write(rootGenome)
    
    if globals.printToConsole:
        #Output newick tree after the ancestors have been added to it
        Phylo.draw(newickTree)
    
    #Need to traverse tree to ouput appropriate content to file
    newickTree.clade.name = '' #Make sure that the output for the root is not output
    traverseNewickTreeAndOutputToFile(newickTree.clade)
        
    #Output the totals for the computation to console and file
    outputTotalsToFile(outputFileName, totalTime)
    
    #TODO compute lineage
    #target = 'NC_014019'
    #print('Computing lineage cost for: %s' % (target))
    #lineageCost = computeLineageCost(newickTree.clade, target, None)
    #if lineageCost != None:
        #print('Successfully found and computed the lineage for: %s' % (target))
    
    if globals.printToConsole:
        #Output Bar graphs of each event
        createBarGraph(globals.deletionSizeCounter, 'Distribution of Deletions')
        createBarGraph(globals.duplicationSizeCounter, 'Distribution of Duplications')
        createBarGraph(globals.inversionSizeDistributionCounter, 'Distribution of Inversions')
        createBarGraph(globals.transpositionSizeDistributionCounter, 'Distribution of Transpositions')
        createBarGraph(globals.invertedTranspositionSizeDistributionCounter, 'Distribution of Inverted Transpositions')
        
    print('Total time (in seconds): %s' % (totalTime))
    print('Ending application...')
    
main()