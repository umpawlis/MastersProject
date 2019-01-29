import time
import os.path
from Bio import Phylo
from strain import Strain

#Parameters that user will pass in
newickFileName = 'Bacillus_Tree.dnd'

#Global variables used in script
strains = [] #Global variable that stores all of the strains in the phylogeny
ancestralCounter = 0 #Counter used to create a unique name for the ancestor

#################################################
########Functions used in this script############
#################################################

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

        if neighborStrain == None:
            print('No neighbor found!')
        else:
            print('The following neighbor will be used during the comparison: %s' % (neighborStrain.getName()))

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

    if neighbor is None and len(currNode.clades) > 0:
        neighbor = getNeighborStrain(currNode.clades[0])

    if neighbor is None and len(currNode.clades) > 1:
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
startTime = time.time()

print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree recursively reconstructing ancestral genomes
result = traverseNewickTree(newickTree.clade, None)

endTime = time.time()
print('Total number of seconds: %s\nProcessing Complete!' % (endTime - startTime))