import time
import os.path
from Bio import Phylo
from strain import Strain
from FileService import createFile

newickFileName = 'Bacillus_Tree.dnd' #Name of newick tree file

#Global variables used in script
strains = [] #Global variable that stores all of the strains in the phylogenetic tree
ancestralCounter = 0 #Counter used to create a unique name for the ancestral node
outputFileName = 'ApplicationOutput.txt'

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
# traverseNewickTree
# Parameters: node - Strain being currently processed, parentNode - direct ancestor of node
# Description: Traverses a provided newick tree in post order traversal
######################################################
def traverseNewickTree(node, parentNode):
    global strains
    global ancestralCounter
    
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
            print('Retrieving strain from strains list: %s' % (foundStrain.getName()))
            return foundStrain
        else:
            newStrain = createStrainFromFile(node)
            if newStrain != None:
                print('Successfully created the following strain from data file: %s' % (newStrain.getName()))
                strains.append(newStrain)
                return newStrain


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