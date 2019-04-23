import time
import os.path
import globals
from Bio import Phylo
from FileService import createFile
from FileService import processSequence
from GlobalAlignmentModule import findOrthologsByGlobalAlignment
from SelfGlobalAlignmentModule import findOrthologsBySelfGlobalAlignment

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
    ancestor = None
    
    print('Computing orthologous events: %s, %s' % (strain1.name, strain2.name))
    orthologousEvents = constructEvents(strain1, strain2)
    
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
    coverageTracker1 = initializeTracker(strain1) #Tracks which operons have been marked off in strain 1
    coverageTracker2 = initializeTracker(strain2) #Tracks which operons have been marked off in strain 2
    
    print('Performing global alignment with: %s, %s' % (strain1.name, strain2.name))
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)
    
    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))
    
    #TODO Add Local Alignment
    
    #Self Global Alignment
    if numRemainingOperons1 > 0:
        duplicationEvents1, lossEvents1, coverageTracker1, strain1 = findOrthologsBySelfGlobalAlignment(strain1, coverageTracker1)
        print('%s, duplicates identified %s and losses identified %s' % (strain1.name, len(duplicationEvents1), len(lossEvents1)))
        if len(lossEvents1) > 0:
            events.extend(lossEvents1)
            
    if numRemainingOperons2 > 0:
        duplicationEvents2, lossEvents2, coverageTracker2, strain2 = findOrthologsBySelfGlobalAlignment(strain2, coverageTracker2)
        print('%s, duplicates identified %s and losses identified %s' % (strain2.name, len(duplicationEvents2), len(lossEvents2)))
        if len(lossEvents2) > 0:
            events.extend(lossEvents2)
            
    #Verify there's no unmarked operons at this point
    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    if numRemainingOperons1 > 0 or numRemainingOperons2 > 0:
        print('Error! There are unmarked operons remaining!')
    
    print(strain1.name)
    print(strain1.codonMismatchDetails)
    print(strain1.substitutionDetails)
    print(strain1.duplicationDetails)
    print(strain1.deletionDetails)
    
    print(strain2.name)
    print(strain2.codonMismatchDetails)
    print(strain2.substitutionDetails)
    print(strain2.duplicationDetails)
    print(strain2.deletionDetails)
    
    return events
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

globals.initialize()
createFile(outputFileName) #Creates file where data will be output

print('Reading newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree recursively reconstructing ancestral genomes
print('Traversing newick tree...')
result = traverseNewickTree(newickTree.clade, None)

endTime = time.time()
totalTime = endTime - startTime
print('Total time (in seconds): %s' % (totalTime))

print('Ending application...')