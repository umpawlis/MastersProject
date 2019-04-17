import time
import os.path
import globals
from Bio import Phylo
from FileService import createFile
from FileService import appendToFile
from FileService import processSequence
from GlobalAlignmentModule2 import findOrthologsByGlobalAlignment

newickFileName = 'Bacillus_Tree.dnd' #Name of newick tree file

#Global variables used in script
strains = [] #Global variable that stores all of the strains in the phylogenetic tree
outputFileName = 'ApplicationOutput.txt'

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
# Description: Identifies all of the orthologous operons between the strains being commpared
######################################################
def constructEvents(strain1, strain2):
    events = []
    coverageTracker1 = {}
    coverageTracker2 = {}
    
    sequence1 = strain1.genomeFragments #Genome fragments of left sibling
    sequence2 = strain2.genomeFragments #Genome fragments of right sibling
    
    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False
    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False

    #Global Alignment operation
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)
    print('Number of orthologous operons identified using Global Alignment %s' % (globalAlignmentCounter))

    numRemainingOperons1 = countRemainingOperons(coverageTracker1)
    numRemainingOperons2 = countRemainingOperons(coverageTracker2)
    print('The number of remaining operons in each respective tracker is: %s, %s' % (numRemainingOperons1, numRemainingOperons2))
    
    return events, strain1, strain2

######################################################
# processStrains
# Parameters:
# Description: Takes two sibling strains and a neighbor and attempts to reconstruct the ancestor
######################################################
def processStrains(strain1, strain2, neighborStrain):
    ancestralSequence = []
    events = []
    
    print('Computing orthologous operons for strains: %s, %s' % (strain1, strain2))
    events, strain1, strain2 = constructEvents(strain1, strain2)

    #Add content to output file
    appendStrainToFile(strain1)
    appendStrainToFile(strain2)
    
    return None

######################################################
# appendStrainToFile
# Parameters:
# Description: Adds content to the output file based on the strain's data
######################################################
def appendStrainToFile(strain):
    appendToFile(outputFileName, "%s\n" % ("Strain:" + strain.name)) #Name of the strain
    appendToFile(outputFileName, "%s\n" % (strain.codonMismatchDetails)) #Codon Mismatches
    appendToFile(outputFileName, "%s\n" % (strain.duplicationInfo)) #Duplication
    appendToFile(outputFileName, "%s\n" % (strain.deletionInfo)) #Deletion
    
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

        ancestor = processStrains(leftSibling, rightSibling, neighborStrain)

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