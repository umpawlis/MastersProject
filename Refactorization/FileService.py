import copy
from Bio import Phylo
from GenomeFragment import GenomeFragment
from BacterialStrain import BacterialStrain

####################################################
######File Service Functions Refactorization########
####################################################

######################################################
# outputStrainDetailsToFile
# Parameters:
# Description: Creates a file where the output will be stored
######################################################
def outputStrainDetailsToFile(fileName, strain,  inversionDetails, transpositionDetails, invertedTransposedDetails):
    #Append all details to file here
    appendToFile(fileName, 'Strain:' + strain.name + '\n')
    appendToFile(fileName, strain.codonMismatchDetails + '\n')
    appendToFile(fileName, strain.substitutionDetails + '\n')
    appendToFile(fileName, strain.duplicationDetails + '\n')
    appendToFile(fileName, strain.deletionDetails + '\n')
    appendToFile(fileName, inversionDetails + '\n')
    appendToFile(fileName, transpositionDetails + '\n')
    appendToFile(fileName, invertedTransposedDetails + '\n')
    
    print(strain.name)
    print(strain.codonMismatchDetails)
    print(strain.substitutionDetails)
    print(strain.duplicationDetails)
    print(strain.deletionDetails)
    print(inversionDetails)
    print(transpositionDetails)
    print(invertedTransposedDetails)

######################################################
# createFile
# Parameters:
# Description: Creates a file where the output will be stored
######################################################
def createFile(fileName, newickTree):
    print('Creating file %s...' % (fileName));    
    file = open(fileName, "w+")
    Phylo.write(newickTree,fileName, 'newick') #Write the tree to the output file
    file.close()
    
######################################################
# appendToFile
# Parameters:
# Description: Appends content to file
######################################################
def appendToFile(fileName, content):
    print('Appending to file %s the following content, %s...' % (fileName, content))
    file = open(fileName, "a")
    file.write(content)
    file.close()

######################################################
# processSequence
# Parameters: sequence - The strain's genome sequence
# Description: Parses the genome sequence into a strain object
######################################################
def processSequence(name, genome):
    index = 0
    genomePosition = 0
    fragmentIndex = 0
    fragments = []
    
    while index < len(genome):
        
        #Case 1: It's either a < o > or a < t >
        if genome[index] == '<':
            startIndex = index
            originOrTerminusPosition = genomePosition
            
            while genome[index] != '>':
                index+=1
            index+=1 #Increment the index to include the >
            genomePosition+=1 #Increment the genome position
            
            originOrTerminus = genome[startIndex:index]
            
            if 't' in originOrTerminus:
                description = 'Terminus'
            else:
                description = 'Origin'
            
            fragment = GenomeFragment(fragmentIndex, originOrTerminus, originOrTerminus.split(','), originOrTerminusPosition, description, False)
            fragments.append(fragment)
            fragmentIndex+=1
        
        #Case 2: It's an operon
        elif (genome[index] == '[') or (genome[index] == '-' and genome[index + 1] == '['):
            startIndex = index
            operonStart = genomePosition
            
            while genome[index] != ']':
                if genome[index] == ',':
                    genomePosition+=1
                index+=1
            index+=1 #Increment index to include the ]
            genomePosition+=1 #Increment to include the last gene when we hit the ]
            
            operon = genome[startIndex:index]
            negativeOrientation = '-' in operon
            genes = []
            
            #Process the operon removing any un-necessary characters and split based on the comma
            data = copy.deepcopy(operon)
            data = data.replace('-', '')
            data = data.replace('[', '')
            data = data.replace(']', '')
            data = data.split(',')
            
            for gene in data:
                genes.append(gene.strip())
            
            if negativeOrientation: #Make positive orientation if negative so all operons are in the same orientation when performing the alignment
                genes.reverse()
            
            fragment = GenomeFragment(fragmentIndex, operon, genes, operonStart, 'Operon', negativeOrientation)
            fragments.append(fragment)
            fragmentIndex+=1
            
        #Case 3: It's a singleton gene
        elif (genome[index] == '-' and genome[index+1].isalnum()) or (genome[index].isalnum()):
            startIndex = index
            singletonStart = genomePosition
            
            while index < len(genome) and genome[index] != ',':
                index+=1
            genomePosition+=1
            
            singleton = genome[startIndex:index]
            negativeOrientation = '-' in singleton
            
            #Remove the negative signed chatacter if it is present
            data = copy.deepcopy(singleton)
            data = data.replace('-', '')
            genes = data.split(',')
            
            fragment = GenomeFragment(fragmentIndex, singleton, genes, singletonStart, 'Singleton', negativeOrientation)
            fragments.append(fragment)
            fragmentIndex+=1
            
        index+=1
        
    return BacterialStrain(name, fragments)