from Bio import Phylo
import os.path
import multiset
import numpy as np
import xlsxwriter

ancestralCounter = 0
threshold = 2
fullAlignmentCounter = 0
extensionCounter = 0

######################################################
# Tracking Event
# Parameters: 
# Description: Stores information about a pair of orthologous operons
######################################################
class TrackingEvent(object):
    trackingEventId = 0
    score = -1
    genome1Name = ""
    genome2Name = ""
    genome1Operon = ""
    genome2Operon = ""
    genome1OperonIndex = -1
    genome2OperonIndex = -1
    ancestralOperon = ""
    technique = ""
    
    def __init__(self, trackingEventId, score, genome1Name, genome2Name, genome1Operon, genome2Operon, genome1OperonIndex, genome2OperonIndex, ancestralOperon, technique):
        self.trackingEventId = trackingEventId
        self.score = score
        self.genome1Name = genome1Name
        self.genome2Name = genome2Name
        self.genome1Operon = genome1Operon
        self.genome2Operon = genome2Operon
        self.genome1OperonIndex = genome1OperonIndex
        self.genome2OperonIndex = genome2OperonIndex
        self.ancestralOperon = ancestralOperon
        self.technique = technique
        
    def printTrackingEvent(self):
        print("{ Tracking Event Id: %d, Alignment Score: %d, Genome 1: %s, Genome 2: %s, Genome 1 Operon: %s, Genome 2 Operon: %s, Genome 1 Operon Index: %d, Genome 2 Operon Index: %d, Ancestral Operon: %s, Technique: %s}"%(self.trackingEventId, self.score, self.genome1Name, self.genome2Name, self.genome1Operon, self.genome2Operon, self.genome1OperonIndex, self.genome2OperonIndex, self.ancestralOperon, self.technique))

    #####################################
    #Getters
    #####################################
    def getTrackingEventId(self):
        return self.trackingEventId
    
    def getScore(self):
        return self.score
    
    def getGenome1Name(self):
        return self.genome1Name
    
    def getGenome2Name(self):
        return self.genome2Name
    
    def getGenome1Operon(self):
        return self.genome1Operon
    
    def getGenome2Operon(self):
        return self.genome2Operon
    
    def getGenome1OperonIndex(self):
        return self.genome1OperonIndex
    
    def getGenome2OperonIndex(self):
        return self.genome2OperonIndex
    
    def getAncestralOperon(self):
        return self.ancestralOperon
    
    def getTechnique(self):
        return self.technique
    
    #####################################
    #Setters
    #####################################
    def setTrackingEventId(self, trackingEventId):
        self.trackingEventId = trackingEventId
        
    def setScore(self, score):
        self.score = score
        
    def setGenome1Name(self, genome1Name):
        self.genome1Name = genome1Name
        
    def setGenome2Name(self, genome2Name):
        self.genome2Name = genome2Name
        
    def setGenome1Operon(self, genome1Operon):
        self.genome1Operon = genome1Operon
    
    def setGenome2Operon(self, genome2Operon):
        self.genome2Operon = genome2Operon
        
    def setGenome1OperonIndex(self, genome1OperonIndex):
        self.genome1OperonIndex = genome1OperonIndex
        
    def setGenome2OperonIndex(self, genome2OperonIndex):
        self.genome2OperonIndex = genome2OperonIndex
        
    def setAncestralOperon(self, ancestralOperon):
        self.ancestralOperon = ancestralOperon
        
    def setTechnique(self, technique):
        self.technique = technique

######################################################
# Strain
# Parameters: name, sequence, descendants
# Description: Stores information about a Strain from the tree
######################################################
class Strain(object):
    name = ""
    sequence = []
    genes = []
    descendants = []
    operonPositions = []
    singletonDict = {}
    
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
        
######################################################
# computeSetDifference
# Parameters: operon1, operon2 - the two operons to compare
# Description: computes the set difference between two operons
# and processes the operon lists
######################################################
def computeSetDifference(operon1, operon2):
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
    
    for op in operon1List:
        set1.add(op.split('_')[0].strip())
        noDuplicatesSet1.add(op.split('_')[0].strip())
    
    for op in operon2List:
        set2.add(op.split('_')[0].strip())
        noDuplicatesSet2.add(op.split('_')[0].strip())
    
    set3 = set1.symmetric_difference(set2)
    noDuplicatesSet3 = noDuplicatesSet1.symmetric_difference(noDuplicatesSet2)
    
    return len(set3), operon1List, operon2List, len(noDuplicatesSet3)
    
######################################################
# sequenceAnalysis
# Parameters: directoryList - List of directories to process
# Description: 
###################################################### 
def sequenceAnalysis(firstOperonList, secondOperonList, strain1, strain2, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    
    print('Starting sequence analysis of {%s, %s}...' % (strain1, strain2))
    
    #remove unwanted content from the operon list
    firstOperonList = removeContent(firstOperonList)
    secondOperonList = removeContent(secondOperonList)
    
    #initialize the matrix to store the global alignment scores
    globalAlignmentMatrix = [[ 0 for x in range(0, len(secondOperonList))] for y in range(0, len(firstOperonList))]
    
    ####################################
    ##Calculations
    ####################################
    for x in range(0, len(firstOperonList)):
        for y in range(0, len(secondOperonList)):
            op1 = firstOperonList[x]
            op2 = secondOperonList[y]
                    
            #check if we need to reverse any of the operons
            reverseOp1 = reverseSequence(op1)
            reverseOp2 = reverseSequence(op2)
            
            #compute the set differences between the two operons
            setDifference, operon1, operon2, numDifferentGenes = computeSetDifference(op1, op2)
            
            #Reverse operons if needed to
            if reverseOp1:
                operon1.reverse()
            
            if reverseOp2:
                operon2.reverse()
                
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
                            scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + 0.5
                    else:
                        scoreMatrix[a][b] = min(scoreMatrix[a-1][b] + 1, scoreMatrix[a][b-1] + 1, scoreMatrix[a-1][b-1] + 1)
                        
            globalAlignmentMatrix[x][y] = scoreMatrix[len(operon1)][len(operon2)]
            
            #comment out if threshold of two is desired, use 1/3 of longest operon as threshold
            threshold = max(len(operon1), len(operon2))
            threshold = threshold//3
            
            if setDifference <= threshold:
                globalAlignmentMatrix[x][y] = str(globalAlignmentMatrix[x][y]) + '*'
                
    ####################################
    ##End of Calculations
    ####################################
    print ('Done analyzing {%s, %s}\n' % (strain1, strain2)) 
    printStrains(strain1, strain2)
    outputResults(strain1, strain2, firstOperonList, secondOperonList, globalAlignmentMatrix)
    outputResultsToExcel(strain1, strain2, firstOperonList, secondOperonList, globalAlignmentMatrix)
    ancestralOperons = findOrthologs(strain1, strain2, firstOperonList, secondOperonList, globalAlignmentMatrix, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)
     
    return ancestralOperons

######################################################
# findOrthologs
# Parameters: strain1, strain2, sequence1, sequence2, resultMatrix
# Description: Scans the matrix to find orthologs
######################################################
def findOrthologs(strain1, strain2, sequence1, sequence2, resultMatrix, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    
    ancestralOperons = []
    
    #Initialize the arrays that will track the coverage of the operons in the two genomes
    coverageTracker1 = {}
    coverageTracker2 = {}
    
    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False
        
    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False
    
    #Use global alignment scores to find orthologs
    coverageTracker1, coverageTracker2, ancestralOperons, numGlobalAlignment, numLocalAlignment = findOrthologsWithGlobalAlignment(strain1, strain2, resultMatrix, coverageTracker1, coverageTracker2, sequence1, sequence2, ancestralOperons, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)
    
    print('#####################################################################')
    print('Statistics for the following strains: %s, %s' %(strain1, strain2))
    print('Number of orthologs found through global alignment: %d' %(numGlobalAlignment))
    print('Number of orthologs found through local alignment: %d' %(numLocalAlignment))
    print('#####################################################################')
          
    ##########################
    # Printer for debugging
    ##########################
    print('Remaining operons from each respective tracker')
    for x in range(0, len(coverageTracker1)):
        if coverageTracker1[x] == False:
            print ('Sequence 1, index: %d, Operon: %s' % (x, sequence1[x]))
            
    for x in range (0, len(coverageTracker2)):
        if coverageTracker2[x] == False:
            print('Sequence 2, index: %d, Operon: %s' % (x, sequence2[x]))
    print('Fin\n')
    
    return ancestralOperons

######################################################
# findOrthologsWithGlobalAlignment
# Parameters: globalAlignmentMatrix, coverageTracker1, coverageTracker2, sequence1, sequence2
# Description: Finds orthologous operons using global alignment
######################################################
def findOrthologsWithGlobalAlignment(genomeName1, genomeName2, globalAlignmentMatrix, coverageTracker1, coverageTracker2, sequence1, sequence2, ancestralOperons, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    
    #Keep track of number of global and local alignments
    globalAlignmentCounter = 0
    localAlignmentCounter = 0
    
    #Tracking Events store information about the ortholog
    trackingId = 0
    trackingEvents = []
    
    #Scan each row
    for i in range(0, len(globalAlignmentMatrix)):
        
        #Track the lowest score
        lowestScore = -1
        
        #Scan each entry in a row
        for j in range(0, len(globalAlignmentMatrix[i])):   
            
            #Check if the entry has an asterisk and if both operons have not been marked off
            if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False):
                currentScore = float(str(globalAlignmentMatrix[i][j]).replace("*", ""))
                
                #If we have not found a score yet or this score is lower then select it
                if lowestScore == -1 or currentScore < lowestScore:
                    lowestScore = currentScore
                    rowIndex = i
                    colIndex = j
                    distance = abs(i - j)
                
                #If we have a score that is the same, select the one closest to the diagonal
                elif lowestScore == currentScore and (abs(i - j) < distance):
                    lowestScore = currentScore
                    rowIndex = i
                    colIndex = j
                    distance = abs(i - j)

        #If we found an ortholog, then mark off both operons
        if lowestScore > -1:
            globalAlignmentCounter+=1
            coverageTracker1[rowIndex] = True
            coverageTracker2[colIndex] = True
            ancestralOperons.append(sequence1[rowIndex])
            
            trackingId += 1
            trackingEvent = TrackingEvent(trackingId, lowestScore, genomeName1, genomeName2, sequence1[rowIndex], sequence2[colIndex], rowIndex, colIndex, "", "2 Genome Global Alignment")
            trackingEvents.append(trackingEvent)
            
            #Used for debugging
            #print('Found an orthologous operon using Global Alignment: (left of matrix) %s, (top of matrix) %s' %(sequence1[rowIndex], sequence2[colIndex]))
            #print('These are the indexes of the orthologous operon from the global alignment: (left of matrix) %d, (top of matrix) %d\n' %(rowIndex, colIndex))
        #Else we failed to find a suitable global alignment so perform a local alignment
        else:
            localAlignmentCounter+=1
            
            op1 = sequence1[i]

            print('\n\n**************************************')
            print('**************************************\n')
            for j in range(0, len(globalAlignmentMatrix[i])):
                if coverageTracker1[i] == False and coverageTracker2[j] == False:
                    op2 = sequence2[j]
                    localAlignment(op1, op2, i, j, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)
            print('\n**************************************')
            print('**************************************\n\n')
    printStats()
    
    if len(trackingEvents) > 0:
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        for i in range(0, len(trackingEvents)):
            trackingEvents[i].printTrackingEvent()
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    
    return coverageTracker1, coverageTracker2, ancestralOperons, globalAlignmentCounter, localAlignmentCounter

def printStrains(strain1, strain2):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nResults of Local Alignment: %s & %s\n" %(strain1, strain2))
    file.write("...\n")
    file.close()

def printAlignments(op1Position, op2Position, operon1, operon2, alignment1, alignment2, message):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("Strain 1 Operon %d: %s\n" %(op1Position, operon1))
    file.write("Strain 2 Operon  %d: %s\n" %(op2Position, operon2))
    file.write("Alignment Result (%s):\n" %message)
    file.write("%s\n" %alignment1)
    file.write("%s\n\n" %alignment2)

    file.close()

def printStats():
    global fullAlignmentCounter
    global extensionCounter
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nTotals:\n")
    file.write("Number of full alignments that occurred: %d\n" %fullAlignmentCounter)
    file.write("Number of extensions performed: %d\n\n" %extensionCounter)

    file.close()
    fullAlignmentCounter = 0
    extensionCounter = 0

######################################################
# localAlignment
# Parameters: operon1, operon2
# Description: performs local alignment on two operons
######################################################
def localAlignment(op1, op2, op1Position, op2Position, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    global fullAlignmentCounter
    print('\nPerforming local alignment..')
    print(op1)
    print(op2)

    #algorithm parameters
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5
    mismatch = -1.0
    gap = -1.0
            
    #check if we need to reverse any of the operons
    reverseOp1 = reverseSequence(op1)
    reverseOp2 = reverseSequence(op2)
    
    #compute the set differences between the two operons
    setDifference, operon1, operon2, numberOfDifferentGenes = computeSetDifference(op1, op2)
    
    leftAdjustment1 = 0
    leftAdjustment2 = 0
    rightAdjustment1 = 0
    rightAdjustment2 = 0
    #Reverse operons if needed to
    if reverseOp1:
        operon1.reverse()
    
    if reverseOp2:
        operon2.reverse()

    #initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))
    maxScore = 0
    maxPosition = (0, 0)

    #perform the Local Alignment
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):

            #check if genes are identical
            if operon1[a-1].split('_')[0].strip() == operon2[b-1].split('_')[0].strip():

                #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                if operon1[a-1].strip() == operon2[b-1].strip():
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + matchWithCodon
                else:
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + matchWithoutCodon
            else:
                scoreMatrix[a][b] = max(scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)

            if scoreMatrix[a][b] > maxScore:
                maxScore = scoreMatrix[a][b]
                maxPosition = (a, b)

    aligned1, aligned2, numGaps, endPosition = traceback(operon1, operon2, scoreMatrix, maxPosition)
    #print(scoreMatrix)

    if len(operon1) <= len(operon2):
        shortestLength = len(operon1)
    else:
        shortestLength = len(operon2)

    if reverseOp1:
        leftAdjustment1 = len(operon1)-1
        leftAdjustment1 -= endPosition[0]-1
        rightAdjustment1 += (len(operon1) - (maxPosition[0]))
    else:
        leftAdjustment1 += endPosition[0]-1
        rightAdjustment1 = len(operon1)-1
        rightAdjustment1 -= (len(operon1) - (maxPosition[0]))
    
    if reverseOp2:
        leftAdjustment2 = len(operon2)-1
        leftAdjustment2 -= endPosition[1]-1
        rightAdjustment2 += (len(operon2) - (maxPosition[1]))
    else:
        leftAdjustment2 += endPosition[1]-1
        rightAdjustment2 = len(operon2)-1
        rightAdjustment2 -= (len(operon2) - (maxPosition[1]))

    if numGaps == 0 and genesStrain1 and genesStrain2:
        if len(aligned1) == shortestLength:
            fullAlignmentCounter += 1
            print("One operon is a SUBSET of the other. Trying extention..")
            extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left and right extension")
        elif len(aligned1) < shortestLength:
            if maxPosition[0] == len(scoreMatrix)-1 and maxPosition[1] == len(scoreMatrix[0])-1:
                print("Trying right extension..")
                extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after right extension")
            elif endPosition[0] == 1 and endPosition[1] == 1:
                print("Trying left extension..")
                extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left extension")
    elif numGaps == 0 and len(aligned1) == shortestLength:
        fullAlignmentCounter += 1
        printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after no extension due to missing gene list(s)")
    print("Final alignment:")
    print(aligned1)
    print(aligned2)

def extendAlignment(direction, operon1, operon2, genesStrain1, genesStrain2, opGenePosition1, opGenePosition2, leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2):
    global extensionCounter
    mismatch = False
    operonRange1 = range(opGenePosition1, opGenePosition1+len(operon1))
    operonRange2 = range(opGenePosition2, opGenePosition2+len(operon2))
    opGenePosition1 += leftAdjustment1
    opGenePosition2 += leftAdjustment2

    while not mismatch:
        if direction == "left":
            if reverseOp1:
                opGenePosition1 += 1
            else:
                opGenePosition1 -= 1

            if reverseOp2:
                opGenePosition2 += 1
            else:
                opGenePosition2 -= 1
        elif direction == "right":
            if reverseOp1:
                opGenePosition1 -= 1
            else:
                opGenePosition1 += 1

            if reverseOp2:
                opGenePosition2 -= 1
            else:
                opGenePosition2 += 1

        #print(str(opGenePosition1) + " " + str(opGenePosition2))
        if (opGenePosition1 in range(0, len(genesStrain1))) and (opGenePosition2 in range(0, len(genesStrain2))):
            if (str(opGenePosition1) in singletonDict1 and str(opGenePosition2) in singletonDict2):
                if (('-' in singletonDict1[str(opGenePosition1)]) != reverseOp1) or (('-' in singletonDict2[str(opGenePosition2)]) != reverseOp2):
                    mismatch = True
            elif (str(opGenePosition1) in singletonDict1 and opGenePosition2 in operonRange2):
                if (('-' in singletonDict1[str(opGenePosition1)]) != reverseOp1):
                    mismatch = True
            elif (str(opGenePosition2) in singletonDict2 and opGenePosition1 in operonRange1):
                if (('-' in singletonDict2[str(opGenePosition2)]) != reverseOp2):
                    mismatch = True
            else:
                mismatch = True

            if not mismatch:
                gene1 = genesStrain1[opGenePosition1]
                gene2  = genesStrain2[opGenePosition2]
                #print(gene1 + " " + gene2)

                if gene1 == gene2:
                    aligned1.insert(0, gene1)
                    aligned2.insert(0, gene2)
                    extensionCounter += 1
                else:
                    mismatch = True
        else:
            mismatch = True

######################################################
# traceback
# Parameters: operon1, operon2, scoreMatrix, startPosition
# Description: treaverses score matrix to determine optimal alignment
######################################################
def traceback(operon1, operon2, scoreMatrix, startPosition):

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = startPosition
    move = nextMove(scoreMatrix, x, y)
    numGaps = 0

    while move != END:
        if move == DIAG:
            aligned_seq1.append(operon1[x - 1])
            aligned_seq2.append(operon2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(operon1[x - 1])
            aligned_seq2.append('-')
            x -= 1
            numGaps += 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(operon2[y - 1])
            y -= 1
            numGaps += 1

        move = nextMove(scoreMatrix, x, y)

    aligned_seq1.append(operon1[x - 1])
    aligned_seq2.append(operon2[y - 1])
    endPosition = (x,y)

    return list(reversed(aligned_seq1)), list(reversed(aligned_seq2)), numGaps, endPosition

######################################################
# nextMove
# Parameters: scoreMatrix, x, y
# Description: determines which direction to traverse score matrix
######################################################
def nextMove(scoreMatrix, x, y):
    diag = scoreMatrix[x - 1][y - 1]
    up   = scoreMatrix[x - 1][y]
    left = scoreMatrix[x][y - 1]

    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')

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
# outputResults
# Parameters: strain1, strain2, sequence1, sequence2, resultMatrix
# Description: outputs the results into a file
######################################################
def outputResults(stain1, strain2, sequence1, sequence2, resultMatrix):
    f = open('out.txt','a')
    
    f.write('%s : %s' % (stain1, sequence1))
    f.write('%s : %s' % (strain2, sequence2))
    f.write('\n')
    f.write('This is the distance matrix for %s %s' % (stain1, strain2))
    f.write('\n')
    
    for row in resultMatrix:
        for val in row:
            f.write(str(val) + ' \t')
        f.write('\n')
    
    f.write('\n')
    f.close()
    
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
# removeContent
# Parameters: operonList - list of operons
# Description: removes unwanted content from the operon list
######################################################
def removeContent(operonList):
    
    if '< o >' in operonList:
        operonList.remove('< o >')
        
    if '< t >' in operonList:
        operonList.remove('< t >')
        
    if '<originCycling>' in operonList:
        operonList.remove('<originCycling>')
        
    return operonList

######################################################
# getOperons
# Parameters: sequence - sequence to get the operons from
# Description: extracts a list of operons from a sequence
###################################################### 
def getOperons(sequence):
    operonList = []
    singletonList = {}
    index = 0
    geneIndex = 0
    operonPositions = []
    
    while index < len(sequence):
        
        #get the origin and termus (remove if not needed)
        if sequence[index] == '<':
            startIndex = index
            
            while sequence[index] != '>':
                index += 1
            
            #increment the index to include the >
            index += 1
            operonList.append(sequence[startIndex:index])
            #decrement the geneIndex to not include <> values
            geneIndex -= 1
            
        elif (sequence[index] == '[') or (sequence[index] == '-' and sequence[index + 1] == '['):
            startIndex = index
            operonPositions.append(geneIndex)
            
            while sequence[index] != ']':
                if sequence[index] == ',':
                    geneIndex += 1
                index += 1
            
            #increment the index to include the ]
            index += 1
            operonList.append(sequence[startIndex:index])

        if sequence[index] == ',':
            geneIndex += 1
            if sequence[index+2] != '[' and sequence[index+2] != '<' and sequence[index+3] != '[':
                index += 1
                startIndex = index

                while sequence[index] != ',' and index < len(sequence)-1:
                    index += 1

                singletonList[str(geneIndex)] = sequence[startIndex:index]
                index -= 1
        index += 1
    #print(sequence)
    #print(operonList)
    #print(operonPositions)
    #print(singletonList)
    return operonList, operonPositions, singletonList

######################################################
# processDistanceFile
# Parameters: fileName - Name of the file to process
# Description: Reads in the file and stores the genes as an array
######################################################
def processDistanceFile(fileName):
    genes1 = []
    genes2 = []
    listTwo = False
    
    file = open(fileName, 'r')
    
    for line in file:   
        #Split the data
        data = line.split()
        #The third position holds the gene
        if 'tRNA' in data[0] or 'rRNA' in data[0]:
            gene = data[2]
            #Remove the tRNA-
            if 'tRNA-' in gene:
                gene = gene.replace("tRNA-", "")
            #add the gene to the genes array
            if listTwo:
                genes2.append(gene)
            else:
                genes1.append(gene)
        elif 'Origin' in data[0]:
            listTwo = True

    genes2.extend(genes1)

    #Close the file
    file.close()
    return genes2
 
######################################################
# post_traversal
# Parameters: node - The node that we want to process
# Description: iterates to left and right nodes, reads 
# in the sequence and takes the symmetric difference between the sequences
###################################################### 
def post_traversal(node):
    
    leftChildStrain = None
    rightChildStrain = None
    currNodeOperons = None
    
    #Check if the clade has children
    if len(node.clades) > 0:
        leftChildStrain = post_traversal(node.clades[0])
        rightChildStrain = post_traversal(node.clades[1])
    
    #Check if the clade has a name, if it does, check if it has a directory for its sequence
    if node.name is not None and len(node.name) > 0:
        print('\nProcessing the node: %s' % node.name)
        
        if os.path.isdir(node.name):
            print('There exists a directory for the node: %s' % node.name)
            
            if os.path.isfile(node.name + '/sequence.rtf'):
                print('Opening the file: %s/sequence.rtf' % node.name)
                
                #Read the sequence in
                fileGeneSequence = open(node.name + '/sequence.rtf', 'r').read()
                
                #Get the operons for this sequence
                currNodeOperons, operonPositions, singletonDict = getOperons(fileGeneSequence)
                strain = Strain(node.name, currNodeOperons, [], operonPositions, singletonDict)                
                
                #Check if we have the distance tRNA & rRNA distance file for this strain
                if os.path.isfile(node.name + '/' + node.name + '_tRNA_and_rRNA_Positions.txt'):
                    print('Opening tRNA rRNA distance file: %s/%s_tRNA_and_rRNA_Positions.txt' % (node.name, node.name))
                    
                    #Read the distance file in
                    allGenes = processDistanceFile(node.name + '/' + node.name + '_tRNA_and_rRNA_Positions.txt')
                    strain.setGenes(allGenes)
                    
                return strain
    
    if leftChildStrain is not None and leftChildStrain.getSequence() is not None and len(leftChildStrain.getSequence()) > 0 and rightChildStrain is not None and rightChildStrain.getSequence() is not None and len(rightChildStrain.getSequence()) > 0:
        
        print('\nThese are the strains that will be compared:')
        leftChildStrain.printStrain()
        rightChildStrain.printStrain()
        
        ancestralOperons = sequenceAnalysis(leftChildStrain.getSequence(), rightChildStrain.getSequence(), leftChildStrain.getName(), rightChildStrain.getName(), leftChildStrain.getGenes(), rightChildStrain.getGenes(), leftChildStrain.getOperonPositions(), rightChildStrain.getOperonPositions(), leftChildStrain.getSingletonDict(), rightChildStrain.getSingletonDict())
        
        global ancestralCounter
        ancestralCounter += 1
        
        node.name = 'Ancestor %d' % (ancestralCounter)

        #TO DO: Properly calculate operon positions for the ancestor
        ancestor = Strain('Ancestor %d' % (ancestralCounter), ancestralOperons, [leftChildStrain.getName(), rightChildStrain.getName()], [], {})
        
        print('This is the resulting ancestor after the comparison:')
        ancestor.printStrain()
        
        return ancestor
        
    #If the left child has a sequence, return it
    elif leftChildStrain is not None and leftChildStrain.getSequence() is not None and len(leftChildStrain.getSequence()) > 0:
        return leftChildStrain
    
    #If the right child has a sequence, return it
    elif rightChildStrain is not None and rightChildStrain.getSequence() is not None and len(rightChildStrain.getSequence()) > 0:
        return rightChildStrain
    
    #If neither has a sequence, return None
    else:
        return None        
        
######################################################
#                       main
###################################################### 
print 'Reading in phylogenetic tree...'
tree = Phylo.read('simpletree.dnd', 'newick')
print 'Done reading in phylogenetic tree'

open('localAlignmentResults.txt', 'w+').close()
#Traverses, computes the ancestral genomes and returns the root node
result = post_traversal(tree.clade)

#Output the root node
if result is not None:
    print('This is the result:')
    result.printStrain()

#Draw tree to the console
Phylo.draw(tree)

print 'End of processing'