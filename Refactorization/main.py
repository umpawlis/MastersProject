import time
import os.path
from Bio import Phylo
from FileService import createFile
from FileService import processSequence

#Application parameters
newickFileName = 'Bacillus_Tree.dnd' #Name of newick tree file

#Global vaiables used by the script
strains = [] #Global variable that stores all of the strains in the phylogenetic tree
outputFileName = 'ApplicationOutput.txt' #Name of output file

################################################################
################# Main Fuctions ################################
################################################################

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

        #ancestor = processStrains(leftSibling, rightSibling, neighborStrain)

        return None
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