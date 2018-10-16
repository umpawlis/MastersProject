from Bio import Phylo
import os.path
import multiset
import numpy as np
import xlsxwriter
import matplotlib.pyplot as plt

newickFileName = 'Anc27_subtree.dnd'
strains = []
ancestralCounter = 0
deletionCost = 1
substitutionCost = 2
codonCost = 0.5

######################################################
# Strain
# Parameters:
# Description: Stores information about a Strain from the newick tree
######################################################
class Strain(object):
    name = ""
    sequence = []
    genes = []
    descendants = []
    operonPositions = []
    singletonDict = {}
    trackingEvents = []
    hasData = False

    #Class constructor
    def __init__(self, name, sequence, descendants, operonPositions, singletonDict):
        self.name = name
        self.sequence = sequence
        self.descendants = descendants
        self.operonPositions = operonPositions
        self.singletonDict = singletonDict

    #Prints the Strain content
    def printStrain(self):
        print("{ Name: %s, Sequence: %s, Descendants: %s }" %(self.name, self.sequence, self.descendants))

    def getName(self):
        return self.name
    def getSequence(self):
        return self.sequence
    def getDescendants(self):
        return self.descendants
    def getGenes(self):
        return self.genes
    def setGenes(self, genes):
        self.genes = genes
    def getOperonPositions(self):
        return self.operonPositions
    def getSingletonDict(self):
        return self.singletonDict
    def getTrackingEvents(self):
        return self.trackingEvents
    def setTrackingEvents(self, trackingEvents):
        self.trackingEvents = trackingEvents
    def setHasData(self, hasData):
        self.hasData = hasData
    def getHasData(self):
        return self.hasData

######################################################
# processFileSequence
# Parameters: sequence - sequence to get the operons from
# Description: extracts a list of operons from a sequence
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
# traverseNewickTree
# Parameters: node - The node that we want to process
# Description: Traverses a provided newick tree and reads in sequences if any (uses post traversal)
######################################################
def traverseNewickTree(node):
    #Global variables
    global strains
    global ancestralCounter
    #Local variables
    leftChildStrain = None
    rightChildStrain = None

    #Check if the clade has children
    if len(node.clades) > 0:
        leftChildStrain = traverseNewickTree(node.clades[0])
        if len(node.clades) > 1:
            rightChildStrain = traverseNewickTree(node.clades[1])

    #Check if the clade has a name, if it does, check if it has a directory for its sequence
    if node.name is not None and len(node.name) > 0:
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
                strain.setHasData(True)
                strains.append(strain)

                return strain

            else:
                print('No sequence file found for node: %s' % node.name)
        else:
            print('No directory found for node: %s' % node.name)

    if leftChildStrain is not None and leftChildStrain.getHasData() and rightChildStrain is not None and rightChildStrain.getHasData():
        print('These are the strains being compared: %s, %s'%(leftChildStrain.getName(), rightChildStrain.getName()))
        ancestralOperons, trackingEvents = processStrains(leftChildStrain, rightChildStrain)

        ancestralCounter += 1
        node.name = 'Ancestor %s' % (ancestralCounter)
        ancestor = Strain('Ancestor %s' % (ancestralCounter), ancestralOperons, [leftChildStrain.getName(), rightChildStrain.getName()], [], {})
        ancestor.setTrackingEvents(trackingEvents)
        ancestor.setHasData(False)
        strains.append(ancestor)

        return ancestor

    #If the left child has a sequence, return it
    elif leftChildStrain is not None and leftChildStrain.getHasData():
        return leftChildStrain

    #If the right child has a sequence, return it
    elif rightChildStrain is not None and rightChildStrain.getHasData():
        return rightChildStrain

    #If neither has a sequence, return None
    else:
        return None

######################################################
# processStrains
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def processStrains(strain1, strain2):
    coverageTracker1 = {}
    coverageTracker2 = {}
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()

    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False

    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False

    #TODO Alignments
    detectOrthologsByGlobalAlignment(strain1, strain2)

    trackerDebugger(coverageTracker1, coverageTracker2, sequence1, sequence2)

    return None, None

######################################################
# processStrains
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def detectOrthologsByGlobalAlignment(strain1, strain2):
    
    #Compute the global alignment matrix
    globalAlignmentMatrix, operonEventMatrix = computeGlobalAlignmentMatrix(strain1, strain2)
    
    
    return None

######################################################
# computeGlobalAlignmentMatrix
# Parameters: directoryList - List of directories to process
# Description:
######################################################
def computeGlobalAlignmentMatrix(strain1, strain2):
    
    firstOperonList = strain1.getSequence()
    secondOperonList = strain2.getSequence()
    strain1Name = strain1.getName()
    strain2Name = strain2.getName()
    
    print('Computing global alignment matrix for: {%s, %s}...' % (strain1Name, strain2Name))

    #initialize the matrix to store the global alignment scores
    globalAlignmentMatrix = [[ 0.0 for x in range(0, len(secondOperonList))] for y in range(0, len(firstOperonList))]
    operonEventMatrix = [[None for x in range(0, len(secondOperonList))] for y in range(0, len(firstOperonList))]

    ####################################
    ##Calculations
    ####################################
    for x in range(0, len(firstOperonList)):
        for y in range(0, len(secondOperonList)):
            op1 = firstOperonList[x]
            op2 = secondOperonList[y]

            #compute the set differences between the two operons
            setDifference, operon1, operon2, numDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

            #Checks if either operons at in the - orientation
            reverseOp1 = reverseSequence(op1)
            reverseOp2 = reverseSequence(op2)

            #Reverse operons if needed to
            if reverseOp1 and reverseOp2 == False:
                operon1.reverse()
            if reverseOp2 and reverseOp1 == False:
                operon2.reverse()

            #Case 1: We have two singleton operons
            if len(operon1) == 1 and len(operon2) == 1:
                #Perfect match
                if operon1 == operon2:
                    globalAlignmentMatrix[x][y] =  str(0) + '*'
                #Mismatch
                else:
                    globalAlignmentMatrix[x][y] = -999

            #Case 2: Only one of them is a singleton operon
            elif (len(operon1) == 1 and len(operon2) > 1) or (len(operon2) == 1 and len(operon1) > 1):
                globalAlignmentMatrix[x][y] = -999

            #Case 3: None of them are singleton operons, perform a global alignment
            elif len(op1) > 1 and len(op2) > 1:
                score, operonEvents = performGlobalAlignment(operon1, operon2)

                globalAlignmentMatrix[x][y] = score
                operonEventMatrix[x][y] = operonEvents

                threshold = max(len(operon1), len(operon2))
                threshold = threshold//3

                if setDifference <= threshold:
                    globalAlignmentMatrix[x][y] = str(globalAlignmentMatrix[x][y]) + '*'

            #Case 4: Some unhandled case
            else:
                print('Case 4: Error, an unhandled case has occured in the sequence analysis')

    ####################################
    ##End of Calculations
    ####################################
    print ('Done computing global alignment matrix for {%s, %s}\n' % (strain1Name, strain2Name))
    outputResultsToExcel(strain1Name, strain2Name, firstOperonList, secondOperonList, globalAlignmentMatrix)
    
    return globalAlignmentMatrix, operonEventMatrix

######################################################
# outputResultsToExcel
# Parameters: strain1, strain2, sequence1, sequence2, resultMatrix
# Description: outputs the results into an excel file
######################################################
def outputResultsToExcel(strain1, strain2, sequence1, sequence2, resultMatrix):

    rowIndex = 1

    #Intialize the workbook
    workbook = xlsxwriter.Workbook('excel_files/Comparison - %s & %s.xlsx' % (strain1, strain2))
    worksheet = workbook.add_worksheet()
    worksheet2 = workbook.add_worksheet()

    #Used to emphasize potentially interesting operons
    cellFormat = workbook.add_format()
    cellFormat.set_pattern(1)
    cellFormat.set_bg_color('lime')

    #Write strain identifiers to the excel file
    worksheet.write(5, 0, strain1)
    worksheet.write(0, 5, strain2)
    worksheet2.write(5, 0, strain1)
    worksheet2.write(0, 5, strain2)

    #Write the operons to the excel file
    for x in range (0, len(sequence1)):
        worksheet.write((x + 2), 1, sequence1[x])
        worksheet2.write((x + 2), 1, sequence1[x])

    for x in range (0, len(sequence2)):
        worksheet.write(1, (x + 2), sequence2[x])
        worksheet2.write(1, (x + 2), sequence2[x])

    #Write the data from the matrix to the excel file
    for row in resultMatrix:

        #Track the excel indexes
        rowIndex += 1
        colIndex = 2

        for value in row:
            if '*' in str(value):
                worksheet.write(rowIndex, colIndex, str(value), cellFormat)
            else:
                worksheet.write(rowIndex, colIndex, str(value))
            colIndex += 1

    #Close the excel file
    workbook.close()

######################################################
# reverseSequence
# Parameters: operon - operon to check
# Description: checks if the operon needs to be reversed
######################################################
def reverseSequence(operon):

    reverseOperon = False

    if '-' in operon:
        reverseOperon = True

    return reverseOperon

######################################################
# formatAndComputeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: computes the set difference between two operons
# and processes the operon lists
######################################################
def formatAndComputeOperonDifferences(operon1, operon2):
    noDuplicatesSet1 = set()
    noDuplicatesSet2 = set()

    set1 = multiset.Multiset();
    set2 = multiset.Multiset();

    operon1 = operon1.replace('-', '')
    operon1 = operon1.replace('[', '')
    operon1 = operon1.replace(']', '')

    operon2 = operon2.replace('-', '')
    operon2 = operon2.replace('[', '')
    operon2 = operon2.replace(']', '')

    operon1List = operon1.split(',')
    operon2List = operon2.split(',')

    noWhiteSpaceOperon1List = []
    noWhiteSpaceOperon2List = []

    for op in operon1List:
        set1.add(op.split('_')[0].strip())
        noDuplicatesSet1.add(op.split('_')[0].strip())
        noWhiteSpaceOperon1List.append(op.strip())

    for op in operon2List:
        set2.add(op.split('_')[0].strip())
        noDuplicatesSet2.add(op.split('_')[0].strip())
        noWhiteSpaceOperon2List.append(op.strip())

    set3 = set1.symmetric_difference(set2)
    noDuplicatesSet3 = noDuplicatesSet1.symmetric_difference(noDuplicatesSet2)

    return len(set3), noWhiteSpaceOperon1List, noWhiteSpaceOperon2List, len(noDuplicatesSet3)

######################################################
# trackerDebugger
# Parameters:
# Description: Prints the content of each respective tracker
######################################################
def trackerDebugger(coverageTracker1, coverageTracker2, sequence1, sequence2):
    print('Remaining operons from each respective tracker:')
    for x in range(0, len(coverageTracker1)):
        if coverageTracker1[x] == False:
            print ('Sequence 1, index: %s, Operon: %s' % (x, sequence1[x]))
    for x in range (0, len(coverageTracker2)):
        if coverageTracker2[x] == False:
            print('Sequence 2, index: %s, Operon: %s' % (x, sequence2[x]))
    print('Finished printing trackers\n')
    
######################################################
# performGlobalAlignment
# Parameters:
# Description: Performs a forward and reversed global alignment and returns the best score
######################################################
def performGlobalAlignment(operon1, operon2):
    global substitutionCost
    global deletionCost
    global codonCost

    #initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))

    for a in range(0, len(operon1)+1):
        scoreMatrix[a][0] = a

    for a in range(0, len(operon2)+1):
        scoreMatrix[0][a] = a

    #perform the Global Alignment
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            #check if genes are identical
            if operon1[a-1].split('_')[0].strip() == operon2[b-1].split('_')[0].strip():
                #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                if operon1[a-1].strip() == operon2[b-1].strip():
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1]
                else:
                    #Solves a special case with a bunch of duplicates with different codons
                    scoreMatrix[a][b] = min(scoreMatrix[a-1][b-1] + codonCost, scoreMatrix[a-1][b] + deletionCost, scoreMatrix[a][b-1] + deletionCost, scoreMatrix[a-1][b-1] + substitutionCost)
            else:
                scoreMatrix[a][b] = min(scoreMatrix[a-1][b] + deletionCost, scoreMatrix[a][b-1] + deletionCost, scoreMatrix[a-1][b-1] + substitutionCost)

    #Compute the number of events that occured between the operons
    operonEvents = globalAlignmentTraceback(scoreMatrix, operon1, operon2)

    return scoreMatrix[len(operon1)][len(operon2)], operonEvents

######################################################
# globalAlignmentTraceback
# Parameters:
# Description: Performs a traceback on a given matrix
######################################################
def globalAlignmentTraceback(matrix, operon1, operon2):
    global substitutionCost
    global deletionCost
    global codonCost

    i = len(operon1)
    j = len(operon2)

    match = 0
    codonMismatch = 0
    mismatch = 0
    substitution = 0
    
    operon1Gaps = []
    operon1Gap = []
    operon1ConsecutiveGap = False #Tracks consecutive gaps

    operon2Gaps = []
    operon2Gap = []
    operon2ConsecutiveGap = False #Tracks consecutive gaps

    #Track the alignment
    alignmentSequence1 = []
    alignmentSequence2 = []

    while i > 0 or j > 0:
        #Case 1: Perfect match
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] and operon1[i-1] == operon2[j-1]:
            match += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + codonCost) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            codonMismatch += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 3: Substitution
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + substitutionCost):
            substitution += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and matrix[i][j] == (matrix[i-1][j] + deletionCost):
            index = i-1
            mismatch += 1
            i -= 1
            operon1ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon2ConsecutiveGap:
                operon2Gap.append(operon1[index])
                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.append(operon2Gap)
                operon2Gap = []
                operon2Gap.append(operon1[index])
                operon2ConsecutiveGap = True

        #Case 5: Mismatch - Gap in operon 1
        else:
            index = j - 1
            mismatch += 1
            j -= 1
            operon2ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon1ConsecutiveGap:
                operon1Gap.append(operon2[index])
                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.append(operon1Gap)
                operon1Gap = []
                operon1Gap.append(operon2[index])
                operon1ConsecutiveGap = True
    
    #Empty any remaining gaps                
    if len(operon1Gap) > 0:
        operon1Gaps.append(operon1Gap)
        operon1Gap = []
    if len(operon2Gap) > 0:
        operon2Gaps.append(operon2Gap)
        operon2Gap = []
        
    #Computes number of unique genes and if the genes are not unique then they are removed from the gap
    operon1Gaps, geneDuplicateSizesInOperon1 = checkForMatchesInAlignment(operon1Gaps, alignmentSequence1)
    operon2Gaps, geneDuplicateSizesInOperon2 = checkForMatchesInAlignment(operon2Gaps, alignmentSequence2)
    
    #operonEvents = OperonEvents(match, codonMismatch, mismatch, substitution, operon1, operon2, matrix, operon1Losses, operon2Losses, operon1Duplications, operon2Duplications)

    #Sets the information we need to perform the sliding window later
    #operonEvents.setAlignedGenesInOperon1(alignmentSequence1)
    #operonEvents.setAlignedGenesInOperon2(alignmentSequence2)
    #operonEvents.setOperon1Gaps(operon1Gaps)
    #operonEvents.setOperon2Gaps(operon2Gaps)
    #operonEvents.setDuplicateSizesOp1(geneDuplicateSizesInOperon1)
    #operonEvents.setDuplicateSizesOp2(geneDuplicateSizesInOperon2)
    
    #return operonEvents
    
######################################################
# checkForMatchesInAlignment
# Parameters:
# Description: Takes an array of gaps and an alignment, then checks if the genes in the gap match any of the genes in the alignment, if they do then genes are popped off the gap array
######################################################
def checkForMatchesInAlignment(arrayOfGaps, alignedGenes):
    geneDuplicateSizes = []
    
    for gap in arrayOfGaps:
        #Initialize Window
        windowSize = len(gap)
        startIndex = 0
        endIndex = len(gap)
        
        #print('Current gap %s' % (gap))
        while windowSize > 1:
            genes = gap[startIndex:endIndex]
            
            #print('Current Window %s' %(genes))
            genesMatched = 0
            
            for x in range(0, len(alignedGenes)):
                if len(genes) > 0 and genes[0] == alignedGenes[x] and genesMatched == 0:
                    for y in range(0, len(genes)):
                        if (x+y) < len(alignedGenes) and genes[y] == alignedGenes[x+y]:
                            genesMatched +=1
                    if genesMatched != len(genes):
                        genesMatched = 0
                        
            if genesMatched != 0 and genesMatched == len(genes):
                #print("Duplicate")
                #updateDuplicationCounter(len(genes))
                geneDuplicateSizes.append(len(genes))
                del gap[startIndex:endIndex]                
                startIndex = endIndex 
            else:
                startIndex+=1
                
            if (startIndex + windowSize) > len(gap):
                #reduce and reset
                windowSize = min(windowSize-1, len(gap))
                startIndex = 0
            endIndex = startIndex + windowSize

    return arrayOfGaps, geneDuplicateSizes

######################################################
#                       main
######################################################
print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree to reconstruct the ancestral genomes
result = traverseNewickTree(newickTree.clade)