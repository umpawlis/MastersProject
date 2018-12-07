from Bio import Phylo
from strain import Strain
from operonEvent import OperonEvents
from trackingEvents import TrackingEvent
import os.path
import multiset
import numpy as np
import xlsxwriter
import matplotlib.pyplot as plt
import copy
import time

#Parameters for script
newickFileName = 'Bacillus_Tree.dnd'
debug = True

#Global Variables
strains = []
ancestralCounter = 0
deletionCost = 1
substitutionCost = 1
codonCost = 0.5
trackingId = 0
duplicateOperonCounter = {}
deletionEventCounter = {}
duplicationEventCounter = {}
yDistanceThreshold = 3

######################################################
# updateDuplicationCounter
# Parameters: duplicationSize
# Description: Updates the duplication counter when a duplication event occurs.
######################################################
def updateDuplicationCounter(duplicationSize):
    global duplicationEventCounter

    if str(duplicationSize) in duplicationEventCounter:
        duplicationEventCounter[str(duplicationSize)] += 1
    else:
        duplicationEventCounter[str(duplicationSize)] = 1

######################################################
# checkOverlap
# Parameters: newRange, rangeList
# Description: Checks if the newly matched genes were already matched within another set.
######################################################
def checkOverlap(newRange, rangeList):
    overlap = False

    for ranges in rangeList:
        if (newRange[0] <= ranges[1]) and (ranges[0] <= newRange[1]):
            overlap = True

    return overlap

######################################################
# geneInSequence
# Parameters: gene, sequence, comparisonSize, opIndex
# Description: Determines if the list of genes appear somewhere else in the operon.
######################################################
def geneInSequence(gene, sequence, comparisonSize, opIndex, alignedRange):
    geneFound = False
    currentOpIndex = 0
    currentIndex = 0
    searchOperon = True
    checkRange = False
    rangeList = []
    rangeList.append(alignedRange)

    while currentOpIndex < len(sequence) and not geneFound:
        checkRange = False
        operon = sequence[currentOpIndex]
        if currentOpIndex != opIndex:
            searchOperon = True
        else:
            if alignedRange[0] == 0 and alignedRange[1] == 0:
                searchOperon = False
            else:
                searchOperon = True
                checkRange = True

        if searchOperon:
            while currentIndex+comparisonSize <= len(operon) and not geneFound:
                checkGene = True
                if checkRange:
                    if (currentIndex < alignedRange[0]) or ((currentIndex+comparisonSize-1) > alignedRange[1]):
                        checkGene = False

                if checkGene:
                    operonGene = operon[currentIndex:currentIndex+comparisonSize]
                    if operonGene == gene:
                        geneFound = True

                currentIndex += 1
            currentIndex = 0
        currentOpIndex += 1

    return geneFound

######################################################
# findUniqueGenes
# Parameters: geneList, sequence
# Description: Tries to find the list of genes in another operon of the genome.
######################################################
def findUniqueGenes(geneList, sequence, opIndex, alignedRange=(0,0)):
    comparisonSize = len(geneList)
    startIndex = 0
    currentIndex = 0
    genesNotFound = True
    endOfList = False
    missedGene = False
    newSet = False
    geneRanges = []
    numGeneMatches = 0
    duplicationSizes = []

    while (comparisonSize >= 2) and genesNotFound:
        #Slide the window across the gene list, one gene at a time
        while (startIndex < comparisonSize) and (startIndex+comparisonSize <= len(geneList)) and genesNotFound:

            currentIndex = startIndex
            #Slide the window across the gene list by comparison size
            while currentIndex < len(geneList):
                if currentIndex+comparisonSize <= len(geneList):
                    gene = geneList[currentIndex:currentIndex+comparisonSize]
                    newSet = True
                # elif not missedGene:
                #     gene = geneList[currentIndex:len(geneList)]
                #     newSet = True

                if newSet:
                    if geneInSequence(gene, sequence, len(gene), opIndex, alignedRange):
                        newRange = (currentIndex, currentIndex+comparisonSize-1)
                        if not checkOverlap(newRange, geneRanges):
                            geneRanges.append(newRange)
                            numGeneMatches += len(gene)
                            duplicationSizes.append(len(gene))
                            # updateDuplicationCounter(len(gene))
                        if numGeneMatches == len(geneList):
                            genesNotFound = False
                    else:
                        missedGene = True

                currentIndex += comparisonSize
                newSet = False

            missedGene = False
            startIndex += 1

        currentIndex = 0
        startIndex = 0
        endOfList = False
        #Decrement the window size
        comparisonSize -= 1

    deletionSizes = []
    inSet = True
    deletionSize = 0

    for index in reversed(range(len(geneList))):
        if not checkOverlap((index,index), geneRanges):
            if not inSet:
                inSet = True
            deletionSize += 1
        else:
            geneList.pop(index)
            inSet = False
            if deletionSize != 0:
                deletionSizes.append(deletionSize)
                deletionSize = 0
    #Special case if all genes in list are losses
    if deletionSize != 0:
        deletionSizes.append(deletionSize)
        deletionSize = 0

    return numGeneMatches, deletionSizes, duplicationSizes

######################################################
# updateDeletionCounter
# Parameters: deletionSize
# Description: Updates the deletion counter when a deletion event occurs.
######################################################
def updateDeletionCounter(deletionSize):
    global deletionEventCounter

    if str(deletionSize) in deletionEventCounter:
        deletionEventCounter[str(deletionSize)] += 1
    else:
        deletionEventCounter[str(deletionSize)] = 1

######################################################
# determineAncestor
# Parameters: op1, op2, startPosition, endPosition, aligned1, aligned2
# Description: Determines which ancestor to pick from local alignment results.
######################################################
def determineAncestor(op1, op2, startPosition, endPosition, aligned1, aligned2, sequence1, sequence2, opIndex1, opIndex2, operonEvents):
    chooseOp1 = True
    ancestor = []
    deletionSizes = None
    duplicationSizes = None

    if len(op1) > len(op2):
        unaligned = getUnaligned(op1, startPosition[0]-1, endPosition[0]-1)
        #numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned2, unaligned)
        numUnique = len(op1)-len(aligned2)
        numDuplicateFound = 0
        alignedRange = (startPosition[0]-1, endPosition[0]-1)

        for geneList in unaligned:
            if geneList:
                numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(geneList, sequence1, opIndex1, alignedRange)
                numUnique -= numUniqueFound
                numDuplicateFound += numUniqueFound

        chooseOp1 = False

        if duplicationSizes:
            for size in duplicationSizes:
                updateDuplicationCounter(size)

        operonEvents.setOperon1GeneDuplicates(numDuplicateFound)
        operonEvents.setOperon2GeneLosses(numUnique)

        if deletionSizes:
            for size in deletionSizes:
                updateDeletionCounter(size)

    else:
        unaligned = getUnaligned(op2, startPosition[1]-1, endPosition[1]-1)
        # numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned1, unaligned)
        numUnique = len(op2)-len(aligned1)
        numDuplicateFound = 0
        alignedRange = (startPosition[1]-1, endPosition[1]-1)

        for geneList in unaligned:
            if geneList:
                numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(geneList, sequence2, opIndex2, alignedRange)
                numUnique -= numUniqueFound
                numDuplicateFound += numUniqueFound

        chooseOp1 = False
        operonEvents.setOperon1GeneLosses(numUnique)

        if not (deletionSizes is None):
            for size in deletionSizes:
                updateDeletionCounter(size)

        if not (duplicationSizes is None):
            for size in duplicationSizes:
                updateDuplicationCounter(size)

        operonEvents.setOperon2GeneDuplicates(numDuplicateFound)

    ancestor = unaligned[0] + aligned1 + unaligned[1]

    print(operonEvents.toStringOperonEvents())

    return ancestor

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
# traceback
# Parameters: operon1, operon2, scoreMatrix, startPosition
# Description: treaverses score matrix to determine optimal alignment
######################################################
def traceback(operon1, operon2, scoreMatrix, startPosition):

    match = 1
    codonMismatch = 0
    mismatch = 0
    substitution = 0
    operon1Losses = 0
    operon2Losses = 0
    operon1Duplications = 0
    operon2Duplications = 0

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = startPosition
    move = nextMove(scoreMatrix, x, y)
    numGaps = 0

    while move != END:
        if move == DIAG:
            if operon1[x-1] == operon2[y-1]:
                match += 1
            elif operon1[x-1].split('_')[0].strip() == operon2[y-1].split('_')[0].strip():
                codonMismatch += 1
            else:
                mismatch += 1

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

    operonEvents = OperonEvents(match, codonMismatch, mismatch, substitution, operon1, operon2, scoreMatrix)

    return list(reversed(aligned_seq1)), list(reversed(aligned_seq2)), numGaps, endPosition, operonEvents

######################################################
# printStrains, printAlignments, printStats
# Parameters:
# Description: Prints a summary of the local alignments.
######################################################
def printStrains(strain1, strain2):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nResults of Local Alignment: %s & %s\n" %(strain1, strain2))
    file.write("...\n")
    file.close()

def printAlignments(op1Position, op2Position, operon1, operon2, alignment1, alignment2, message):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("Strain 1 Operon %s: %s\n" %(op1Position, operon1))
    file.write("Strain 2 Operon  %s: %s\n" %(op2Position, operon2))
    file.write("Alignment Result (%s):\n" %message)
    file.write("%s\n" %alignment1)
    file.write("%s\n\n" %alignment2)

    file.close()

def printStats():
    global fullAlignmentCounter
    global extensionCounter
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nTotals:\n")
    file.write("Number of full alignments that occurred: %s\n" %fullAlignmentCounter)
    file.write("Number of extensions performed: %s\n\n" %extensionCounter)

    file.close()
    fullAlignmentCounter = 0
    extensionCounter = 0

######################################################
# extendAlignment
# Parameters: operon1, operon2
# Description: Tries to extend the alignment of two operons by using their gene lists.
######################################################
def extendAlignment(direction, operon1, operon2, genesStrain1, genesStrain2, opGenePosition1, opGenePosition2, leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2):
    # global extensionCounter
    mismatch = False
    operonRange1 = range(opGenePosition1, opGenePosition1+len(operon1))
    operonRange2 = range(opGenePosition2, opGenePosition2+len(operon2))
    opGenePosition1 += leftAdjustment1
    opGenePosition2 += leftAdjustment2

    extensionScore = 0
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5

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

        if (opGenePosition1 in range(0, len(genesStrain1))) and (opGenePosition2 in range(0, len(genesStrain2))):
            #Check if the gene is a singlton or part of an operon
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

                if gene1.split('_')[0].strip() == gene2.split('_')[0].strip():
                    #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                    if gene1.strip() == gene2.strip():
                        extensionScore += matchWithCodon
                    else:
                        extensionScore += matchWithoutCodon

                    aligned1.insert(0, gene1)
                    aligned2.insert(0, gene2)
                    # extensionCounter += 1
                else:
                    mismatch = True
        else:
            mismatch = True
    return extensionScore

######################################################
# getUnaligned
# Parameters: operon, startPosition, endPosition
# Description: Finds genes in the operon that are not part of the alignment
######################################################
def getUnaligned(operon, startPosition, endPosition):
    unaligned = []
    beforeList = []
    afterList = []

    for i in range(len(operon)):
        if i < startPosition:
            beforeList.append(operon[i])
        elif i > endPosition:
            afterList.append(operon[i])

    unaligned.append(beforeList)
    unaligned.append(afterList)

    return unaligned

######################################################
# localAlignment
# Parameters:
# Description: Performs local alignment on two operons
######################################################
def localAlignment(op1, op2, op1Position, op2Position, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    #print('\nPerforming local alignment..')
    #print(op1)
    #print(op2)

    #Algorithm parameters
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5
    mismatch = -1.0
    gap = -1.0

    #Check if we need to reverse any of the operons
    reverseOp1 = reverseSequence(op1)
    reverseOp2 = reverseSequence(op2)

    #Compute the set differences between the two operons
    setDifference, operon1, operon2, numberOfDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

    leftAdjustment1 = 0
    leftAdjustment2 = 0
    rightAdjustment1 = 0
    rightAdjustment2 = 0
    #Reverse operons if needed to
    if reverseOp1:
        operon1.reverse()
    if reverseOp2:
        operon2.reverse()

    #Initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))
    maxScore = 0
    extensionScore = 0
    returningScore = -1
    maxPosition = (0, 0)

    #Perform the Local Alignment
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

    #Trace back score matrix to get alignment
    aligned1, aligned2, numGaps, endPosition, operonEvents = traceback(operon1, operon2, scoreMatrix, maxPosition)

    if len(operon1) <= len(operon2):
        shortestLength = len(operon1)
    else:
        shortestLength = len(operon2)

    #Adjust gene position we are looking at for extension based on the orienation of the operon
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

    #Only extend operons when the current alignment has no gaps
    if numGaps == 0 and genesStrain1 and genesStrain2:
        if len(aligned1) == shortestLength:
            #print("One operon is a SUBSET of the other. Trying extention..")
            extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            returningScore = maxScore + extensionScore
            extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            returningScore += extensionScore
            printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left and right extension")
        elif len(aligned1) < shortestLength:
            if maxPosition[0] == len(scoreMatrix)-1 or maxPosition[1] == len(scoreMatrix[0])-1:
                #print("Trying right extension..")
                extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after right extension")
                    returningScore = maxScore + extensionScore
            elif endPosition[0] == 1 or endPosition[1] == 1:
                #print("Trying left extension..")
                extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left extension")
                    returningScore = maxScore + extensionScore
    #No extension but still qualifies because it has no gaps in alignment
    elif numGaps == 0 and len(aligned1) == shortestLength:
        printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after no extension due to missing gene list(s)")
        returningScore = maxScore
    #print("Final alignment:")
    #print(aligned1)
    #print(aligned2)

    return returningScore, operon1, operon2, maxPosition, endPosition, aligned1, aligned2, operonEvents

######################################################
# formatAllOperons
# Parameters: sequence
# Description: Formats all operons to correct orientation.
######################################################
def formatAllOperons(sequence):
    sequenceList = []
    operonIndexConversions = []
    operonIndex = 0

    for operon1 in sequence:
        reverseOp = reverseSequence(operon1)

        operon1 = operon1.replace('-', '')
        operon1 = operon1.replace('[', '')
        operon1 = operon1.replace(']', '')

        operon1List = operon1.split(',')

        if len(operon1List) > 1:
            noWhiteSpaceOperon1List = []

            for op in operon1List:
                noWhiteSpaceOperon1List.append(op.strip())

            if reverseOp:
                noWhiteSpaceOperon1List.reverse()
            sequenceList.append(noWhiteSpaceOperon1List)
            operonIndexConversions.append(operonIndex)
            operonIndex += 1
        else:
            operonIndexConversions.append(-1)

    return sequenceList, operonIndexConversions

######################################################
# detectOrthologsByLocalAlignment
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def detectOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2):
    global trackingId

    #Finding optimal orthologs using local alignment
    minValue = min(len(coverageTracker1), len(coverageTracker2))
    trackingEvents = []
    localAlignmentCounter = 0
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()
    genomeName1 = strain1.getName()
    genomeName2 = strain2.getName()
    formattedSequence1, operonIndexes1 = formatAllOperons(strain1.getSequence())
    formattedSequence2, operonIndexes2 = formatAllOperons(strain2.getSequence())
    genesStrain1 = strain1.getGenes()
    genesStrain2 = strain2.getGenes()
    operonPositionList1 = strain1.getOperonPositions()
    operonPositionList2 = strain2.getOperonPositions()
    singletonDict1 = strain1.getSingletonDict()
    singletonDict2 = strain2.getSingletonDict()

    #Scan the matrix x times to find optimal local alignments everytime an entire matrix is scanned
    for x in range(0, minValue):
        highestScore = -1
        distance = 50   #arbitrary large number
        chosenOperon1 = []
        chosenOperon2 = []
        chosenStart = (0,0)
        chosenEnd = (0,0)
        chosenAligned1 = []
        chosenAligned2 = []

        for i in range(0, len(coverageTracker1)):
            #Check if operon was not covered and not a singleton
            if coverageTracker1[i] == False and len(sequence1[i].split(',')) > 1:
                for j in range(0, len(coverageTracker2)):
                    if coverageTracker1[i] == False and coverageTracker2[j] == False and len(sequence2[j].split(',')) > 1:
                        op1 = sequence1[i]
                        op2 = sequence2[j]

                        score, formattedOperon1, formattedOperon2, endPosition, startPosition, aligned1, aligned2, operonEvents = localAlignment(op1, op2, i, j, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)

                        if (score > highestScore) or (score == highestScore and (abs(i - j)) < distance):
                            highestScore = score
                            rowIndex = i
                            colIndex = j
                            distance = abs(i - j)

                            chosenOperon1 = formattedOperon1
                            chosenOperon2 = formattedOperon2
                            chosenStart = startPosition
                            chosenEnd = endPosition
                            chosenAligned1 = aligned1
                            chosenAligned2 = aligned2
                            chosenOpEvents = operonEvents

        #After scanning the whole matrix if we found a best score, then store it
        if highestScore > -1:
            print('\n******** Local Alignment ***********')
            print('**************************************\n')

            localAlignmentCounter+=1
            coverageTracker1[rowIndex] = True
            coverageTracker2[colIndex] = True

            ancestralOperon = determineAncestor(chosenOperon1, chosenOperon2, chosenStart, chosenEnd, chosenAligned1, chosenAligned2, formattedSequence1, formattedSequence2, operonIndexes1[rowIndex], operonIndexes2[colIndex], chosenOpEvents)

            trackingId += 1
            trackingEvent = TrackingEvent(trackingId, highestScore, genomeName1, genomeName2, sequence1[rowIndex], sequence2[colIndex], rowIndex, colIndex, ancestralOperon, "Local Alignment")
            trackingEvent.printTrackingEvent()
            trackingEvents.append(trackingEvent)

            print('\n**************************************')
            print('**************************************\n\n')

    return trackingEvents, coverageTracker1, coverageTracker2, localAlignmentCounter

######################################################
# findMaxValueInMatrix
# Parameters:
# Description: Finds the maximum value in the global alignment matrix
######################################################
def findMaxValueInMatrix(globalAlignmentMatrix):

    maxValue = -1
    for i in range(0, len(globalAlignmentMatrix)):
        for j in range(0, len(globalAlignmentMatrix[i])):
            if ('*' in str(globalAlignmentMatrix[i][j])):
                currentValue = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                if currentValue > maxValue:
                    maxValue = currentValue
    return maxValue

######################################################
# scanGlobalAlignmentMatrixForOrthologs
# Parameters:
# Description: Scans matrix and identifies orthologous operons
######################################################
def scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, operonEventMatrix, coverageTracker1, coverageTracker2, strain1, strain2):
    global trackingId
    global codonCost

    genomeName1 = strain1.getName()
    genomeName2 = strain2.getName()
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()
    maxValue = findMaxValueInMatrix(globalAlignmentMatrix)
    currentScoreSelected = 0
    globalAlignmentCounter = 0
    trackingEvents = []

    #Keep iterating util we find all the optimal scores (Finding orthologs using global alignment)
    while currentScoreSelected <= maxValue:
        #Prioritize the selection of operons with the same sign
        for i in range(0, len(globalAlignmentMatrix)):
            for j in range(0, len(globalAlignmentMatrix[i])):
                #Check if this is a * score, if both operons have not been marked off and if both are the same sign
                if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False) and (('-' in sequence1[i] and '-' in sequence2[j]) or ('-' not in sequence1[i] and '-' not in sequence2[j])):
                    score = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                    #Check if the score matches the scores we're currently looking for
                    if score == currentScoreSelected:
                        #We found an ortholog in the global alignment matrix
                        print('\n##### Global Alignment #####')
                        trackingId += 1
                        globalAlignmentCounter+=1
                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        if score == 0:
                            #We found a perfect match, doesn't matter which operon we pick
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, sequence1[i], "2 Genome Global Alignment")
                        else:
                            #We found orthologs that are not a perfect match and need to be resolved
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, '', "2 Genome Global Alignment")
                        #Add the event to the tracking events list
                        trackingEvent.setOperonEvents(operonEventMatrix[i][j])
                        trackingEvents.append(trackingEvent)
                        #TODO
                        trackingEvent.printTrackingEvent()
                        print('###################################\n')

        #Select the remaining operons with the optimal score
        for i in range(0, len(globalAlignmentMatrix)):
            for j in range(0, len(globalAlignmentMatrix[i])):
                #Check if this is a * score and if both operons have not been marked off
                if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False):
                    score = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                    #Check if the score matches the scores we're currently looking for
                    if score == currentScoreSelected:
                        #We found an ortholog in the global alignment matrix
                        print('\n##### Global Alignment #####')
                        trackingId += 1
                        globalAlignmentCounter+=1
                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        if score == 0:
                            #We found a perfect match, doesn't matter which operon we pick
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, sequence1[i], "2 Genome Global Alignment")
                        else:
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, '', "2 Genome Global Alignment")
                        #Add the event to the tracking events list
                        trackingEvent.setOperonEvents(operonEventMatrix[i][j])
                        trackingEvents.append(trackingEvent)
                        trackingEvent.printTrackingEvent()
                        print('###################################\n')
        currentScoreSelected += codonCost

    return trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter

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
# checkForMatchesInAlignment
# Parameters:
# Description: Takes an array of gaps and an alignment, then checks if the genes in the gap match any of the genes in the alignment, if they do then genes are popped off the gap array
# returns an array of duplicate sizes
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

    #Tracks where the extra genes are from
    gap1Indexes = []
    gap2Indexes = []

    #Tracks the substitution
    substitutionDictionary = {}
    operon1IndexToAlignment2Index = {}

    while i > 0 or j > 0:
        #Case 1: Perfect match
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] and operon1[i-1] == operon2[j-1]:
            match += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            operon1IndexToAlignment2Index[str(i-1)] = len(alignmentSequence2)
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + codonCost) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            codonMismatch += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            operon1IndexToAlignment2Index[str(i-1)] = len(alignmentSequence2)
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 3: Substitution
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + substitutionCost):
            substitution += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            operon1IndexToAlignment2Index[str(i-1)] = len(alignmentSequence2)
            substitutionDictionary[str(i-1)] = len(alignmentSequence1)
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and matrix[i][j] == (matrix[i-1][j] + deletionCost):
            index = i-1
            mismatch += 1
            i -= 1
            operon1ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon2ConsecutiveGap:
                operon2Gap.insert(0, operon1[index])
                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.insert(0, operon2Gap)
                operon2Gap = []
                operon2Gap.insert(0, operon1[index])
                gap2Indexes.insert(0, len(alignmentSequence2))
                operon2ConsecutiveGap = True

        #Case 5: Mismatch - Gap in operon 1
        else:
            index = j - 1
            mismatch += 1
            j -= 1
            operon2ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon1ConsecutiveGap:
                operon1Gap.insert(0, operon2[index])
                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.insert(0, operon1Gap)
                operon1Gap = []
                operon1Gap.insert(0, operon2[index])
                gap1Indexes.insert(0, len(alignmentSequence1))
                operon1ConsecutiveGap = True

    #Empty any remaining gaps
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1Gap = []
    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2Gap = []

    #The indexes values need to be flipped b/c right now they're oriented from right to left
    if len(gap1Indexes) > 0:
        for x in range(0, len(gap1Indexes)):
            gap1Indexes[x] = len(alignmentSequence1) - gap1Indexes[x]
    if len(gap2Indexes) > 0:
        for x in range(0, len(gap2Indexes)):
            gap2Indexes[x] = len(alignmentSequence2) - gap2Indexes[x]
    if len(substitutionDictionary) > 0:
        for key, value in substitutionDictionary.items():
            substitutionDictionary[key] = len(alignmentSequence1) - value
    if len(operon1IndexToAlignment2Index) > 0:
        for key, value in operon1IndexToAlignment2Index.items():
            operon1IndexToAlignment2Index[key] = len(alignmentSequence2) - value

    #Need to swap the gap lists since the gaps refer to extra genes
    temp = operon1Gaps
    operon1Gaps = operon2Gaps
    operon2Gaps = temp

    temp = gap1Indexes
    gap1Indexes = gap2Indexes
    gap2Indexes = temp

    #Used for debugging
#    print('These are the operons being compared: %s, %s' %(operon1, operon2))
#    print('This is the resulting alignment: %s, %s' %(alignmentSequence1, alignmentSequence2))
#    print('These are the extra genes for operon 1: %s' %(operon1Gaps))
#    print('These are the indexes for extra genes in operon 1: %s' %(gap1Indexes))
#    print('These are the extra genes for operon 2: %s' %(operon2Gaps))
#    print('These are the indexes for extra genes in operon 2: %s' %(gap2Indexes))
#    print('This is the substitution dictionary: %s' %(substitutionDictionary))
#    print('This is the operon 1 to operon 2 mapper: %s' %(operon1IndexToAlignment2Index))

    operonEvents = OperonEvents(match, codonMismatch, mismatch, substitution, operon1, operon2, matrix)

    #Sets the information we need to perform the sliding window later
    operonEvents.setAlignedGenesInOperon1(alignmentSequence1)
    operonEvents.setAlignedGenesInOperon2(alignmentSequence2)
    operonEvents.setOperon1Gaps(operon1Gaps)
    operonEvents.setOperon2Gaps(operon2Gaps)
    operonEvents.setOperon1GapIndexes(gap1Indexes)
    operonEvents.setOperon2GapIndexes(gap2Indexes)
    operonEvents.setSubstitutionDict(substitutionDictionary)
    operonEvents.setOperon1IndexToAlignment2Index(operon1IndexToAlignment2Index)

    return operonEvents

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
# detectOrthologsByGlobalAlignment
# Parameters:
# Description: Calls function to construct the matrix and then calls function to scan matrix
######################################################
def detectOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2):
    globalAlignmentMatrix, operonEventMatrix = computeGlobalAlignmentMatrix(strain1, strain2)
    trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter = scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, operonEventMatrix, coverageTracker1, coverageTracker2, strain1, strain2)

    return trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter

######################################################
# constructTrackingEvents
# Parameters:
# Description: Constructs the tracking events between two provided strains
######################################################
def constructTrackingEvents(strain1, strain2, printStats):
    coverageTracker1 = {}
    coverageTracker2 = {}
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()

    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False

    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False

    #Global Alignment
    trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter = detectOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)

    #Local Alignment
    localAlignmentTrackingEvents, coverageTracker1, coverageTracker2, localAlignmentCounter = detectOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
    if len(localAlignmentTrackingEvents) > 0:
        trackingEvents.extend(localAlignmentTrackingEvents)

    #Handle singleton genes
    singletonDuplicatedG1, singletonLostG1, singletonTrackingEventsG1, coverageTracker1 = detectOrthologousSingletonGenes(coverageTracker1, strain1, strain1.getTrackingEvents())
    singletonDuplicatedG2, singletonLostG2, singletonTrackingEventsG2, coverageTracker2 = detectOrthologousSingletonGenes(coverageTracker2, strain2, strain2.getTrackingEvents())
    if len(singletonTrackingEventsG1) > 0:
        trackingEvents.extend(singletonTrackingEventsG1)
    if len(singletonTrackingEventsG2) > 0:
        trackingEvents.extend(singletonTrackingEventsG2)

    #Handle remaining operons
    operonDuplicateG1, operonLossG1, operonTrackingEventsG1, coverageTracker1 = detectDuplicateOperons(coverageTracker1, strain1)
    operonDuplicateG2, operonLossG2, operonTrackingEventsG2, coverageTracker2 = detectDuplicateOperons(coverageTracker2, strain2)
    if len(operonTrackingEventsG1) > 0:
        trackingEvents.extend(operonTrackingEventsG1)
    if len(operonTrackingEventsG2) > 0:
        trackingEvents.extend(operonTrackingEventsG2)

    trackerDebugger(coverageTracker1, coverageTracker2, sequence1, sequence2)

    if printStats:
        print('#' * 70)
        print('Statistics for the following strains: %s, %s' %(strain1.getName(), strain2.getName()))
        print('Total number of operons and singletons for %s: %s' %(strain1.getName(), len(coverageTracker1)))
        print('Total number of operons and singletons for %s: %s' %(strain2.getName(), len(coverageTracker2)))
        print('Number of orthologs found through global alignment: %s' %(globalAlignmentCounter))
        print('Number of orthologs found through local alignment: %s' %(localAlignmentCounter))
        print('Number of singletons identified as duplicates in %s: %s' %(strain1.getName(), singletonDuplicatedG1))
        print('Number of singletons identified as duplicates in %s: %s' %(strain2.getName(), singletonDuplicatedG2))
        print('Number of singletons lost in %s: %s' %(strain1.getName(), singletonLostG1))
        print('Number of singletons lost in %s: %s' %(strain2.getName(), singletonLostG1))
        print('Number of operons lost in %s: %s' %(strain1.getName(), operonLossG1))
        print('Number of operons lost in %s: %s' %(strain2.getName(), operonLossG1))
        print('Number of operons identified as duplicates in %s: %s' %(strain1.getName(), operonDuplicateG1))
        print('Number of operons identified as duplicates in %s: %s' %(strain2.getName(), operonDuplicateG2))
        print('#' * 70)

    return trackingEvents

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
# duplicateAlignment
# Parameters:
# Description: Performs a global alignment on genome itself to find best source of duplicate
######################################################
def duplicateAlignment(g1OperonIndex, g1Operon, g1Sequence, genomeName1):
    print('Duplicate Alignment')

    optimalScore = -1
    duplicateEvent = None

    for x in range(0, len(g1Sequence)):
        #Ignore global alignment on operon itself
        if x != g1OperonIndex:
            #Compute the set differences
            setDifference, operon1, operon2, numDifferentGenes = formatAndComputeOperonDifferences(g1Operon, g1Sequence[x])

            #Checks if either operons at in the - orientation
            reverseOp1 = reverseSequence(g1Operon)
            reverseOp2 = reverseSequence(g1Sequence[x])

            #Reverse operons if needed to
            if reverseOp1 and reverseOp2 == False:
                operon1.reverse()
            if reverseOp2 and reverseOp1 == False:
                operon2.reverse()

            #Threshold to check if the sequences are worth comparing
            if setDifference <= (max(len(operon1), len(operon2))//3):
                lowestScore, operonEvents = performGlobalAlignment(operon1, operon2)

                if optimalScore == -1 or (lowestScore < optimalScore):
                    optimalScore = lowestScore
                    duplicateEvent = TrackingEvent(0, optimalScore, genomeName1, genomeName1, g1Operon, g1Sequence[x], g1OperonIndex, x, "", "Duplicate Alignment")
                    duplicateEvent.setOperonEvents(operonEvents)
                    distance = abs(g1OperonIndex - x)
                elif lowestScore == optimalScore and (abs(g1OperonIndex - x) < distance):
                    optimalScore = lowestScore
                    duplicateEvent = TrackingEvent(0, optimalScore, genomeName1, genomeName1, g1Operon, g1Sequence[x], g1OperonIndex, x, "", "Duplicate Alignment")
                    duplicateEvent.setOperonEvents(operonEvents)
                    distance = abs(g1OperonIndex - x)

    #check if we found a duplicate event
    if duplicateEvent:
        incrementDuplicateTracker(operon1)

    return duplicateEvent

######################################################
# incrementDuplicateTracker
# Parameters:
# Description: Increments the tracker according by using the operon length as the key
######################################################
def incrementDuplicateTracker(operon):
    key = str(len(operon)) #Figures out which key we need to increment

    if key in duplicateOperonCounter:
        duplicateOperonCounter[str(key)] = duplicateOperonCounter[str(key)] + 1
    else:
        duplicateOperonCounter[str(key)] = 1

######################################################
# detectDuplicateOperons
# Parameters:
# Description:
######################################################
def detectDuplicateOperons(coverageTracker, strain):
    global trackingId

    sequence = strain.getSequence()
    genomeName = strain.getName()
    operonDuplicate = 0
    operonLoss = 0
    numInvertedDuplicates = 0
    trackingEvents = []

    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False and len(sequence[i].split(',')) > 1:
            print('\n&&&&&&&&&& Duplicate Alignment &&&&&&&&&&&&&&&&')
            print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')

            duplicateEvent = duplicateAlignment(i, sequence[i], sequence, genomeName)
            coverageTracker[i] = True

            #checks if duplicate event is null
            if duplicateEvent:
                #Not null, we found a partner for this operon
                print('Found a matching duplicate operon, therefore not adding to ancestor!')
                duplicateEvent.printTrackingEvent()

                #Increment counter to indicate we successfully found a match
                operonDuplicate += 1

                #Check if inverted duplicates
                if ('-' in duplicateEvent.getGenome1Operon() and '-' not in duplicateEvent.getGenome2Operon()) or ('-' not in duplicateEvent.getGenome1Operon() and '-' in duplicateEvent.getGenome2Operon()):
                    numInvertedDuplicates += 1

            else:
                print('No duplicate ortholog found for operon: %s, therefore it will be added to ancestor as it is a loss' % (sequence[i]))
                trackingId += 1
                trackingEvent = TrackingEvent(trackingId, 0, genomeName, '', sequence[i], '', i, -1, sequence[i], "Duplicate Alignment (No match found)")
                trackingEvent = trackLossEvents(trackingEvent, strain.getTrackingEvents())

                #Indicates operon is a loss
                operonLoss += 1

                #decides whether to add the event or not
                if len(trackingEvent.getLostEventIds()) >= 2:
                    print('Removing this operon because it was lost two times in a row')
                else:
                    trackingEvents.append(trackingEvent)
                trackingEvent.printTrackingEvent()
            print('\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
            print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n')

    return operonDuplicate, operonLoss, trackingEvents, coverageTracker

######################################################
# trackLossEvents
# Parameters:
# Description: Appends the previous tracking loss Ids to the current nodes tracking Id
######################################################
def trackLossEvents(currTrackingEvent, previousTrackingEvents):
    if len(previousTrackingEvents) == 0:
        #If there are no tracking events, then this is a leaf we're dealing with
        lostEventIds = currTrackingEvent.getLostEventIds()
        lostEventIds.append(currTrackingEvent.getTrackingEventId())
        currTrackingEvent.setLostEventIds(lostEventIds)
    else:
        #If there are tracking events, then this is not a leaf. Two cases to consider, it was lost in the previous tracking events or not
        for i in range(0, len(previousTrackingEvents)):
            if previousTrackingEvents[i].getAncestralOperon().strip() == currTrackingEvent.getAncestralOperon().strip():
                #We found the tracking event associated with the operon from the previous ancestral node
                lostEventIds = previousTrackingEvents[i].getLostEventIds()
                lostEventIds.append(currTrackingEvent.getTrackingEventId())
                currTrackingEvent.setLostEventIds(lostEventIds)
                #If there was a loss in the previous ancestor for the same operon, remove the operon from both the previous and current ancestor
                #if len(lostEventIds) >= 2:
                    #previousTrackingEvents.pop(i)

    return currTrackingEvent

######################################################
# resolveSingleton
# Parameters:
# Description: Finds original copy of singleton
######################################################
def resolveSingleton(sequence, singletonIndex, coverageTracker):
    addToAncestor = False
    sourceIndex = -1
    #Check if an exact gene exists in an operon
    for o in range(0, len(sequence)):
        #Don't compare to itself
        if o != singletonIndex:
            #Get a list of operon genes
            setDifference, singletonGene, operonGenes, numDifferentGenes = formatAndComputeOperonDifferences(sequence[singletonIndex], sequence[o])
            #If a match found and no other match found
            if (singletonGene[0] in operonGenes) and sourceIndex == -1:
                sourceIndex = o
                distance = abs(o - singletonIndex)
            #If a match found with smaller distance
            elif (singletonGene[0] in operonGenes) and distance > abs(o - singletonIndex):
                sourceIndex = o
                distance = abs(o - singletonIndex)

    coverageTracker[singletonIndex] = True
    if sourceIndex != -1:
        incrementDuplicateTracker(singletonGene)
        print('\n##### Singleton Source Found!! (NOT ADDED TO ANCESTOR) #####')
        print('The singleton gene %s was found in operon %s, index: %s' % (sequence[singletonIndex].strip(), sequence[sourceIndex].strip(), sourceIndex))
        print('###################################\n')
    else:
        addToAncestor = True
        print('\n##### No Singleton Source Found!!(ADDED TO ANCESTOR) #####')
        print('Could not find source for the singleton gene: %s' % (sequence[singletonIndex].strip()))
        print('###################################\n')

    return addToAncestor, sourceIndex, coverageTracker

######################################################
# detectOrthologousSingletonGenes
# Parameters:
# Description:
######################################################
def detectOrthologousSingletonGenes(coverageTracker, strain, descendantsTrackingEvents):
    global trackingId
    sequence = strain.getSequence()
    genomeName = strain.getName()
    singletonLost = 0
    singletonDuplicated = 0
    trackingEvents = []

    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False and len(sequence[i].split(',')) == 1:
            addToAncestor, matchIndex, coverageTracker = resolveSingleton(sequence, i, coverageTracker)

            if addToAncestor:
                #If no match found then it's a loss so add to ancestor
                trackingId += 1
                trackingEvent = TrackingEvent(trackingId, 0, genomeName, '', sequence[i], '', i, -1, sequence[i], "Singleton Alignment")
                trackingEvent = trackLossEvents(trackingEvent, descendantsTrackingEvents)

                #Indicates a loss
                singletonLost += 1
                #decides whether to add the event or not
                if len(trackingEvent.getLostEventIds()) >= 2:
                    print('Removing this singleton because it was lost two times in a row')
                else:
                    trackingEvents.append(trackingEvent)
                trackingEvent.printTrackingEvent()
            else:
                #Increment Counter to indicate success in finding source
                singletonDuplicated += 1
    return singletonDuplicated, singletonLost, trackingEvents, coverageTracker

######################################################
# processStrains
# Parameters:
# Description: Takes two related strains and a close neighbor and constructs the events for both comparisons
######################################################
def processStrains(strain1, strain2, neighborStrain):

    ancestralSequence = []

    neighborTrackingEvents = None
    if not (neighborStrain is None):
        print('Constructing the tracking events for neighboring strain: %s' % (neighborStrain.getName()))
        neighborTrackingEvents = constructTrackingEvents(strain1, neighborStrain, False)

    print('Constructing tracking events for cherry pair: %s, %s' %(strain1.getName(), strain2.getName()))
    trackingEvents = constructTrackingEvents(strain1, strain2, True)

    if len(trackingEvents) > 0:
        trackingEvents = reconstructAncestralOperon(trackingEvents, strain1, strain2, neighborTrackingEvents)
        if len(neighborTrackingEvents):
            createDotPlot(neighborTrackingEvents, strain1, neighborStrain)
        createDotPlot(trackingEvents, strain1, strain2)
        CFR, TFR, IR, ITR, LO = reconstructAncestralOperonSequence(trackingEvents)
        NCFR, NTFR, NIR, NITR, NLO = reconstructAncestralOperonSequence(neighborTrackingEvents)

        for trackingEvent in trackingEvents:
            stringAncestralOperon = formatAncestralOperontoString(trackingEvent.getAncestralOperon())
            ancestralSequence.append(stringAncestralOperon)
            trackingEvent.setAncestralOperon(stringAncestralOperon)

        #if len(neighborTrackingEvents) > 0:
            #Construct Alignment
            #alignedTrackingEvents = constructAlignment(CFR, TFR, IR, ITR, LO, NCFR, NTFR, NIR, NITR, NLO)
            #Create sequence Array
        #updateGlobalTrackers(trackingEvents, strain1.getSequence(), strain2.getSequence())

    return ancestralSequence, trackingEvents

######################################################
# formatAncestralOperontoString
# Parameters:
# Description: formats the operon into a string
######################################################
def formatAncestralOperontoString(ancestralOperon):
    stringOperon = ""

    if type(ancestralOperon) == str:
        stringOperon = '[' + ancestralOperon + ']'
    elif type(ancestralOperon) == list:
        for x in range(0, len(ancestralOperon)):
            if x == 0:
                stringOperon += '[' + ancestralOperon[x] + ', '
            elif x == len(ancestralOperon) -1:
                stringOperon += ancestralOperon[x] + ']'
            else:
                stringOperon += ancestralOperon[x] + ', '
    else:
        print('Unsupported type detected!')
    return stringOperon

######################################################
# reconstructAncestralOperon
# Parameters:
# Description: Constructs the ancestral operon by determining wheth the extra genes are losses or duplicates
######################################################
def reconstructAncestralOperon(trackingEvents, strain1, strain2, neighborTrackingEvents):
    print('Reconstructing ancestral operons')
    for trackingEvent in trackingEvents:
        if trackingEvent.getTechnique() == '2 Genome Global Alignment' and trackingEvent.getOperonEvents() != None:
            currentOperonEvent = trackingEvent.getOperonEvents()

            if trackingEvent.getScore() == 0:
                #Pefect matches, doesn't matter which one we pick
                trackingEvent.setAncestralOperon(currentOperonEvent.getAlignedGenesInOperon1())
            else:
                #Need to check whether these gaps are losses or duplicates
                operon1Gaps = currentOperonEvent.getOperon1Gaps()
                operon1GapIndexes = currentOperonEvent.getOperon1GapIndexes()
                operon2Gaps = currentOperonEvent.getOperon2Gaps()
                operon2GapIndexes = currentOperonEvent.getOperon2GapIndexes()

                #Used for debugging
#                print('These are the operons being compared: %s, %s' %(trackingEvent.getGenome1Operon(), trackingEvent.getGenome2Operon()))
#                print('This is the resulting alignment: %s, %s' %(currentOperonEvent.getAlignedGenesInOperon1(), currentOperonEvent.getAlignedGenesInOperon2()))
#                print('These are the extra genes for operon 1: %s' %(operon1Gaps))
#                print('These are the indexes for extra genes in operon 1: %s' %(operon1GapIndexes))
#                print('These are the extra genes for operon 2: %s' %(operon2Gaps))
#                print('These are the indexes for extra genes in operon 2: %s' %(operon2GapIndexes))

                #Checks if these extra genes are duplicates by checking if they exist within the alignment and removes them if they do
                operon1Gaps, duplicateSizesWithinAlignment1 = checkForMatchesInAlignment(operon1Gaps, currentOperonEvent.getAlignedGenesInOperon1())
                operon2Gaps, duplicateSizesWithinAlignment2 = checkForMatchesInAlignment(operon2Gaps, currentOperonEvent.getAlignedGenesInOperon2())
                ancestralOperon = currentOperonEvent.getAlignedGenesInOperon1()
                #Decide which ancestral operon to pick
#                subsDict = currentOperonEvent.getSubstitutionDict()
#                if trackingEvent.getScore() > 0 and len(neighborTrackingEvents) > 0 and len(subsDict) > 0:
#                    #Get the neighbors mapper
#                    currentNeighborTrackingEvent = None
#                    for neighborTrackingEvent in neighborTrackingEvents:
#                        if neighborTrackingEvent.getGenome1OperonIndex() == trackingEvent.getGenome1OperonIndex():
#                            currentNeighborTrackingEvent = neighborTrackingEvent
#                            break
#
#                    if currentNeighborTrackingEvent != None:
#                        indexMapper = currentNeighborTrackingEvent.getOperon1IndexToAlignment2Index()
#
#                        ancestralOperon = currentOperonEvent.getAlignedGenesInOperon1()
#
#                        for key, value in subsDict.items():
#                            if indexMapper[key] != None and currentNeighborTrackingEvent.getOperonEvents() != None:
#                                neighborOperonEvent = currentNeighborTrackingEvent.getOperonEvents()
#
#                                neighborAlignment = neighborOperonEvent.getAlignedGenesInOperon2()
#                                neighborGene  = neighborAlignment[indexMapper[key]]
#
#                                operonAlignment = currentOperonEvent.getAlignedGenesInOperon1()
#                                operonGene = operonAlignment[subsDict[key]]
#
#                                #Switch genes if the genes don't match
#                                if operonGene != neighborGene:
#                                    siblingAlignment = currentOperonEvent.getAlignedGenesInOperon2()
#                                    ancestralOperon[subsDict[key]] = siblingAlignment[subsDict[key]]
#                else:
#                    ancestralOperon = currentOperonEvent.getAlignedGenesInOperon1()

                formattedSequence1, operon1SequenceConversion = formatAllOperons(strain1.getSequence())
                formattedSequence2, operon2SequenceConversion = formatAllOperons(strain2.getSequence())
                #print(formattedSequence1)

                i = len(operon1Gaps) - 1
                j = len(operon2Gaps) - 1
                #Testing by inserting more genes
#                if len(operon1Gaps) > 0:
#                    operon1Gaps[0].insert(0, 'Ala_GCA')
#                    operon1Gaps[0].insert(0, 'Gly_GGC')

                while (i > 0) or (j > 0):

                    #Get the biggest index
                    print(i)
                    print(j)
                    if i < 0:
                        #Ran out of elements in Gaps 1
#                        print('Gap being processed: %s' % (operon2Gaps[j]))
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(operon2Gaps[j], formattedSequence2, operon2SequenceConversion[trackingEvent.getGenome2OperonIndex()])
#                        print('Gap being processed: %s' % (operon2Gaps[j]))
#                        print('Number of unique genes found: %s' %(numUniqueFound))
#                        print('Number of deletion genes found: %s' %(deletionSizes))
#                        print('Number of duplicate genes found: %s' %(duplicationSizes))
                        if len(operon2Gaps[j]) > 0:
                            #Insert gap into operon
                            operon2Gaps[j].reverse()
                            for gene in operon2Gaps[j]:
                                ancestralOperon.insert(operon2GapIndexes[j], gene)
                        j = j - 1
                    elif j < 0:
                        #Ran out of elements in Gaps 2
#                        print('Gap being processed: %s' % (operon1Gaps[i]))
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(operon1Gaps[i], formattedSequence1, operon1SequenceConversion[trackingEvent.getGenome1OperonIndex()])
#                        print('Gap being processed: %s' % (operon1Gaps[i]))
#                        print('Number of unique genes found: %s' %(numUniqueFound))
#                        print('Number of deletion genes found: %s' %(deletionSizes))
#                        print('Number of duplicate genes found: %s' %(duplicationSizes))
                        if len(operon1Gaps[i]) > 0:
                            #Insert gap into operon
                            operon1Gaps[i].reverse()
                            for gene in operon1Gaps[i]:
                                ancestralOperon.insert(operon1GapIndexes[i], gene)
                        i = i - 1
                    elif operon1GapIndexes[i] > operon2GapIndexes[j]:
                        #Operon 1 index is bigger
#                        print('Gap being processed: %s' % (operon1Gaps[i]))
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(operon1Gaps[i], formattedSequence1, operon1SequenceConversion[trackingEvent.getGenome1OperonIndex()])
#                        print('Gap being processed: %s' % (operon1Gaps[i]))
#                        print('Number of unique genes found: %s' %(numUniqueFound))
#                        print('Number of deletion genes found: %s' %(deletionSizes))
#                        print('Number of duplicate genes found: %s' %(duplicationSizes))
                        if len(operon1Gaps[i]) > 0:
                            #Insert gap into operon
                            operon1Gaps[i].reverse()
                            for gene in operon1Gaps[i]:
                                ancestralOperon.insert(operon1GapIndexes[i], gene)
                        i = i - 1
                    elif operon1GapIndexes[i] < operon2GapIndexes[j]:
                        #Operon 2 index is bigger
#                        print('Gap being processed: %s' % (operon2Gaps[j]))
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(operon2Gaps[j], formattedSequence2, operon2SequenceConversion[trackingEvent.getGenome2OperonIndex()])
#                        print('Gap being processed: %s' % (operon2Gaps[j]))
#                        print('Number of unique genes found: %s' %(numUniqueFound))
#                        print('Number of deletion genes found: %s' %(deletionSizes))
#                        print('Number of duplicate genes found: %s' %(duplicationSizes))
                        if len(operon2Gaps[j]) > 0:
                            #Insert gap into operon
                            operon2Gaps[j].reverse()
                            for gene in operon2Gaps[j]:
                                ancestralOperon.insert(operon2GapIndexes[j], gene)
                        j = j - 1
                #Set ancestral operon
                trackingEvent.setAncestralOperon(ancestralOperon)
#                print('This is the resulting ancestral operon: %s' % (ancestralOperon))
#                print('\n\n')
#                print('These are the extra genes remaining for operon 1: %s' %(operon1Gaps))
#                print('These are the extra genes remaining for operon 2: %s' %(operon2Gaps))
#                print('These are the duplicate sizes operon 1: %s' %(duplicateSizesWithinAlignment1))
#                print('These are the duplicate sizes operon 2: %s\n\n' %(duplicateSizesWithinAlignment2))

    return trackingEvents

######################################################
# reconstructAncestralOperonSequence
# Parameters:
# Description:
######################################################
def reconstructAncestralOperonSequence(trackingEvents):
    global yDistanceThreshold

    arrayOfConservedForwardRegions = []
    arrayOfTransposedForwardRegions = []
    arrayOfInvertedRegions = []
    arrayOfInvertedTranspositionRegions = []
    arrayOfLostOperons = []

    #Sort Tracking Events on x-coords
    trackingEventsCopy = copy.deepcopy(trackingEvents)
    trackingEventsCopy.sort(key=lambda x: x.genome1OperonIndex, reverse=False)

    while len(trackingEventsCopy) > 0:
        currentTrackingEvent = trackingEventsCopy.pop(0) #Get first item in list

        #Handle lost operons
        if currentTrackingEvent.getGenome1OperonIndex() == -1 or currentTrackingEvent.getGenome2OperonIndex() == -1:
            arrayOfLostOperons.append(currentTrackingEvent)
        else:
            consecutiveRegion = []
            consecutiveRegion.append(currentTrackingEvent)

            foundNeighbor = True
            yIncreaseCount = 0
            yDecreaseCount = 0
            above = False
            below = False

            minMainDiagonalDistance = abs(currentTrackingEvent.getGenome1OperonIndex() - currentTrackingEvent.getGenome2OperonIndex())
            if (currentTrackingEvent.getGenome1OperonIndex() - currentTrackingEvent.getGenome2OperonIndex()) < 0:
                above = True
            if (currentTrackingEvent.getGenome1OperonIndex() - currentTrackingEvent.getGenome2OperonIndex()) > 0:
                below = True

            while foundNeighbor:
                foundNeighbor = False;
                if len(trackingEventsCopy) > 0:
                    previousPoint = consecutiveRegion[len(consecutiveRegion) - 1]
                    currentPoint = trackingEventsCopy[0]
                    distance = abs(previousPoint.getGenome2OperonIndex() - currentPoint.getGenome2OperonIndex())

                    if distance < yDistanceThreshold:
                        #Consecutive point
                        consecutiveRegion.append(trackingEventsCopy.pop(0))
                        foundNeighbor = True

                        currentMainDiagonalDistance = abs(currentPoint.getGenome1OperonIndex() - currentPoint.getGenome2OperonIndex())
                        if currentMainDiagonalDistance < minMainDiagonalDistance:
                            minMainDiagonalDistance = currentMainDiagonalDistance

                        #Checks if main diagonal is crossed
                        if (currentPoint.getGenome1OperonIndex() - currentPoint.getGenome2OperonIndex()) < 0:
                            above = True
                        if (currentPoint.getGenome1OperonIndex() - currentPoint.getGenome2OperonIndex()) > 0:
                            below = True

                        if previousPoint.getGenome2OperonIndex() < currentPoint.getGenome2OperonIndex():
                            yIncreaseCount += 1
                        else:
                            yDecreaseCount +=1

            #Store forward regions
            if yIncreaseCount > yDecreaseCount or (len(consecutiveRegion) == 1):
                if minMainDiagonalDistance < 3:
                    arrayOfConservedForwardRegions.append(consecutiveRegion)
                else:
                    arrayOfTransposedForwardRegions.append(consecutiveRegion)
            #Store Inverted regions
            else:
                if above == True and below == True:
                    arrayOfInvertedRegions.append(consecutiveRegion)
                else:
                    arrayOfInvertedTranspositionRegions.append(consecutiveRegion)
        #print('x-axis: %s, y-axis: %s\n' %(currentTrackingEvent.getGenome1OperonIndex(), currentTrackingEvent.getGenome2OperonIndex()))
    #end while

    print('Stats:')
    print('Total number of tracking events: %s' % (len(trackingEvents)))

    print('Total number of forward conserved regions: %s' % len(arrayOfConservedForwardRegions))
    print('Total number of forward transposed regions: %s' % len(arrayOfTransposedForwardRegions))

    print('Total number of inverted regions: %s' % len(arrayOfInvertedRegions))
    print('Total number of inverted transposed regions: %s' % len(arrayOfInvertedTranspositionRegions))

    print('Total number of lost operons: %s' % len(arrayOfLostOperons))

    for region in arrayOfConservedForwardRegions:
        print('Forward Conserved Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].getGenome1OperonIndex(), region[x].getGenome2OperonIndex()))
    for region in arrayOfTransposedForwardRegions:
        print('Forward Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].getGenome1OperonIndex(), region[x].getGenome2OperonIndex()))

    for region in arrayOfInvertedRegions:
        print('Inverted Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].getGenome1OperonIndex(), region[x].getGenome2OperonIndex()))
    for region in arrayOfInvertedTranspositionRegions:
        print('Inverted Transposed Region')
        for x in range(0, len(region)):
            print('%s, %s' %(region[x].getGenome1OperonIndex(), region[x].getGenome2OperonIndex()))

    for operon in arrayOfLostOperons:
        print('Lost operon')
        print('%s, %s' %(operon.getGenome1OperonIndex(), operon.getGenome2OperonIndex()))

    #ancestralOperonSequence = assembleSequence(arrayOfConservedForwardRegions, arrayOfTransposedForwardRegions, arrayOfInvertedRegions, arrayOfInvertedTranspositionRegions, arrayOfLostOperons)

    return arrayOfConservedForwardRegions, arrayOfTransposedForwardRegions, arrayOfInvertedRegions, arrayOfInvertedTranspositionRegions, arrayOfLostOperons

######################################################
# updateGlobalTrackers
# Parameters:
# Description:
######################################################
def updateGlobalTrackers(trackingEvents, sequence1, sequence2):

    formattedSequence1, operon1SequenceConversion = formatAllOperons(sequence1)
    formattedSequence2, operon2SequenceConversion = formatAllOperons(sequence2)

    for i in range(0, len(trackingEvents)):
        if trackingEvents[i].getTechnique() == '2 Genome Global Alignment':
            opEvents = trackingEvents[i].getOperonEvents()
            #Gets duplicate counts based on similar genes found within alignment
            if opEvents and opEvents.getDuplicateSizesOp1() and len(opEvents.getDuplicateSizesOp1()) > 0:
                for size in opEvents.getDuplicateSizesOp1():
                    updateDuplicationCounter(size)
            if opEvents and opEvents.getDuplicateSizesOp2() and len(opEvents.getDuplicateSizesOp2()) > 0:
                for size in opEvents.getDuplicateSizesOp2():
                    updateDuplicationCounter(size)

            #Computes losses and duplications for each set of gaps
            if opEvents and opEvents.getOperon1Gaps() and len(opEvents.getOperon1Gaps()) > 0:
                gaps = opEvents.getOperon1Gaps()
                totalLosses = 0
                for gap in gaps:
                    if len(gap) > 0:
                        #Compute duplications and losses along with their sizes
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(gap, formattedSequence1, operon2SequenceConversion[trackingEvents[i].getGenome2OperonIndex()])
                        totalLosses += numUniqueFound
                        #Update duplication counter
                        if len(duplicationSizes) > 0:
                            for size in duplicationSizes:
                                updateDuplicationCounter(size)
                        #Update loss counter
                        if len(deletionSizes) > 0:
                            for size in deletionSizes:
                                updateDeletionCounter(size)
                #end for
                opEvents.setLossesDueToSlidingWindowMethodOperon1(totalLosses)

            if opEvents and opEvents.getOperon2Gaps() and len(opEvents.getOperon2Gaps()) > 0:
                gaps = opEvents.getOperon2Gaps()
                totalLosses = 0
                for gap in gaps:
                    if len(gap) > 0:
                        #Uniques tells us number of losses
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(gap, formattedSequence2, operon1SequenceConversion[trackingEvents[i].getGenome1OperonIndex()])
                        totalLosses += numUniqueFound
                        #Update duplication counter
                        if len(duplicationSizes) > 0:
                            for size in duplicationSizes:
                                updateDuplicationCounter(size)
                        #Update loss counter
                        if len(deletionSizes) > 0:
                            for size in deletionSizes:
                                updateDeletionCounter(size)
                #end for
                opEvents.setLossesDueToSlidingWindowMethodOperon2(totalLosses)

######################################################
# createDotPlot
# Parameters:
# Description:
######################################################
def createDotPlot(trackingEvents, strain1, strain2):

    #Stores all of the coordinates
    x_coord = []
    y_coord = []

    #The green ones represent operons with no differences
    green_x_coord = []
    green_y_coord = []

    #Yellow ones represent scores between 1 and 2
    yellow_x_coord = []
    yellow_y_coord = []

    #Red ones represent scores between 3 and above
    orange_x_coord = []
    orange_y_coord = []

    #Blue ones represent a local alignment
    red_x_coord = []
    red_y_coord = []

    print("x" * 70)
    for i in range(0, len(trackingEvents)):
        if trackingEvents[i].getTechnique() == '2 Genome Global Alignment' or trackingEvents[i].getTechnique() == 'Local Alignment':
            #Assign the coords to the appropriate array
            if trackingEvents[i].getTechnique() == 'Local Alignment':
                red_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                red_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            elif trackingEvents[i].getScore() == 0:
                green_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                green_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            elif trackingEvents[i].getScore() == 1 or trackingEvents[i].getScore() == 2:
                yellow_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                yellow_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            else:
                orange_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                orange_y_coord.append(trackingEvents[i].getGenome2OperonIndex())

            #Get all coordinates into a single array
            x_coord.append(trackingEvents[i].getGenome1OperonIndex())
            y_coord.append(trackingEvents[i].getGenome2OperonIndex())

            #print('x-axis: %s, y-axis: %s' %(trackingEvents[i].getGenome1OperonIndex(), trackingEvents[i].getGenome2OperonIndex()))

    #If we have any coordinates to plot, display them
    if len(green_x_coord) > 0 or len(yellow_x_coord) > 0 or len(orange_x_coord) > 0 or len(red_x_coord) > 0:
        f = plt.figure()
        plt.title("Orthologous Operon Mapping")
        plt.plot(green_x_coord, green_y_coord, 'o', color = 'green')
        plt.plot( yellow_x_coord, yellow_y_coord, 'o', color = 'yellow')
        plt.plot(orange_x_coord, orange_y_coord, 'o', color = 'orange')
        plt.plot(red_x_coord, red_y_coord, 'o', color = 'red')
        plt.axis([0, len(trackingEvents)+5, 0, len(trackingEvents)+5])
        plt.ylabel('Operon Position in %s' % (strain1.getName()))
        plt.xlabel('Operon Position in %s' % (strain2.getName()))
        plt.show()
        f.savefig("%s %s.pdf" %(strain1.getName(), strain2.getName()), bbox_inches='tight')
    else:
        print('No plot to display!')
    print("x" * 70)

######################################################
# findNeighboringStrain
# Parameters:
# Description: Finds the first available strain with data
######################################################
def findNeighboringStrain(currNode):
    global strains
    neighbor = None

    if neighbor is None and currNode.name is not None and len(currNode.name) > 0:
        if currNode.is_terminal:
            neighbor = createStrainFromFile(currNode)
        elif 'Ancestor' in currNode.name:
            filteredList = filter(lambda x: x.name == currNode.name, strains)
            neighbor = next(filteredList)
    if neighbor is None and len(currNode.clades) > 0:
        neighbor = findNeighboringStrain(currNode.clades[0])
    if neighbor is None and len(currNode.clades) > 1:
        neighbor = findNeighboringStrain(currNode.clades[1])

    return neighbor

######################################################
# processFileSequence
# Parameters: sequence - sequence to get the operons from
# Description: extracts a list of operons from a sequence file
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
# createStrainFromFile
# Parameters:
# Description: creates strain based on the node provided if sequence data is available
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
# Parameters: node - The node that we want to process
# Description: Traverses a provided newick tree and reads in sequences if any (uses post traversal)
######################################################
def traverseNewickTree(node, parentNode):

    #Global variables
    global strains
    global ancestralCounter

    #Local variables
    leftChildStrain = None
    rightChildStrain = None

    #Check if the clade has children
    if len(node.clades) > 0:
        leftChildStrain = traverseNewickTree(node.clades[0], node)
        if len(node.clades) > 1:
            rightChildStrain = traverseNewickTree(node.clades[1], node)

    #Check if the clade has a name, if it does, create Strain object for it
    if not(node.name == None) and len(node.name) > 0:
        print('%s' % (node.name))
        newStrain = createStrainFromFile(node)
        if not(newStrain == None):
            strains.append(newStrain)
            return newStrain

    if not(leftChildStrain == None) and len(leftChildStrain.getSequence()) > 0 and not(rightChildStrain == None) and len(rightChildStrain.getSequence()) > 0:
        print('These are the strains being compared: %s, %s'%(leftChildStrain.getName(), rightChildStrain.getName()))

        neighborStrain = None
        if parentNode != None:
            node.name = 'Processing'
            #Determine which side we need to traverse
            if len(parentNode.clades) > 0 and parentNode.clades[0].name != "Processing":
                neighborStrain = findNeighboringStrain(parentNode.clades[0])
            elif len(parentNode.clades) > 1 and parentNode.clades[1].name != "Processing":
                neighborStrain = findNeighboringStrain(parentNode.clades[1])

        if debug:
            if neighborStrain == None:
                print('No neighbor found!')
            else:
                print('Found the following neighbor: %s' %(neighborStrain.getName()))

        ancestralOperons, trackingEvents = processStrains(leftChildStrain, rightChildStrain, neighborStrain)

        ancestralCounter += 1
        node.name = 'Ancestor %s' % (ancestralCounter)
        ancestor = Strain('Ancestor %s' % (ancestralCounter), ancestralOperons, [leftChildStrain.getName(), rightChildStrain.getName()], [], {})
        ancestor.setTrackingEvents(trackingEvents)
        strains.append(ancestor)

        return ancestor

    #If the left child has a sequence, return it
    elif not(leftChildStrain == None) and len(leftChildStrain.getSequence()) > 0:
        return leftChildStrain

    #If the right child has a sequence, return it
    elif not(rightChildStrain == None) and len(rightChildStrain.getSequence()) > 0:
        return rightChildStrain

    #If neither has a sequence, return None
    else:
        return None

######################################################
#                       main
######################################################
start = time.time()
print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree to reconstruct the ancestral genomes
result = traverseNewickTree(newickTree.clade, None)

end = time.time()
print(end - start)
print('Done')