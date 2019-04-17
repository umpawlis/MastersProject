import time
from Bio import Phylo
from FileService import createFile


outputFileName = 'ApplicationOutput.txt' #Name of output file
newickFileName = 'Bacillus_Tree.dnd' #Name of newick tree file

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
#result = traverseNewickTree(newickTree.clade, None)

endTime = time.time()
totalTime = endTime - startTime
print('Total time (in seconds): %s' % (totalTime))

print('Ending application...')