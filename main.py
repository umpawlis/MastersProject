import time
import os.path
from Bio import Phylo
from strain import Strain
from GlobalAlignmentModule import findOrthologsByGlobalAlignment
from LocalAlignmentModule import findOrthologsByLocalAlignment
import globals

#Parameters that user will pass in
newickFileName = 'Bacillus_Tree.dnd'

#Global variables used in script
strains = [] #Global variable that stores all of the strains in the phylogeny
ancestralCounter = 0 #Counter used to create a unique name for the ancestor
deletionCost = 1
substitutionCost = 1
codonCost = 0.5


#################################################
########Functions used in this script############
#################################################

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




    return ancestralSequence, events

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

    if numRemainingOperons1 > 0 and numRemainingOperons2 > 0:
        localAlignmentEvents, coverageTracker1, coverageTracker2, localAlignmentCounter = findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
        print('Number of orthologous operons identified using Local Alignment %s' % (localAlignmentCounter))

        numRemainingOperons1 = countRemainingOperons(coverageTracker1)
        numRemainingOperons2 = countRemainingOperons(coverageTracker2)
        print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))

        if len(localAlignmentEvents) > 0:
            events.extend(localAlignmentEvents)

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

    if not(leftSibling == None) and len(leftSibling.getSequence()) > 0 and not(rightSibling == None) and len(rightSibling.getSequence()) > 0:
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

        #TODO: Add the ancestral code
        ancestralOperons, events = processStrains(leftSibling, rightSibling, neighborStrain)

        return None

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
        #Operon
        if (sequence[index] == '[') or (sequence[index] == '-' and sequence[index + 1] == '['):
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
startTime = time.time()

print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree recursively reconstructing ancestral genomes
result = traverseNewickTree(newickTree.clade, None)

endTime = time.time()
print('Total number of seconds: %s\nProcessing Complete!' % (endTime - startTime))