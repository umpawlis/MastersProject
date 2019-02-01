from SequenceService import formatAllOperons
from SequenceService import reverseSequence
from SequenceService import formatAndComputeOperonDifferences
import numpy as np

################################
###Local Alignment Functions####
################################

######################################################
# detectOrthologsByLocalAlignment
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2):
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
            
            #TODO: update this part!
            trackingId += 1
            trackingEvent = TrackingEvent(trackingId, highestScore, genomeName1, genomeName2, sequence1[rowIndex], sequence2[colIndex], rowIndex, colIndex, ancestralOperon, "Local Alignment")
            trackingEvent.printTrackingEvent()
            trackingEvents.append(trackingEvent)

            print('\n**************************************')
            print('**************************************\n\n')

    return trackingEvents, coverageTracker1, coverageTracker2, localAlignmentCounter

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
    operon1NegativeOrientation = reverseOp1
    operon2NegativeOrientation = reverseOp2
    operon1Reversed = False
    operon2Reversed = False
    
    #Compute the set differences between the two operons
    setDifference, operon1, operon2, numberOfDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

    leftAdjustment1 = 0
    leftAdjustment2 = 0
    rightAdjustment1 = 0
    rightAdjustment2 = 0
    #Reverse operons if needed to
    if reverseOp1:
        operon1.reverse()
        operon1Reversed = True
    if reverseOp2:
        operon2.reverse()
        operon2Reversed = True

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
    
    #Store information about the orientation of the operons
    operonEvents.setOperon1NegativeOrientation(operon1NegativeOrientation)
    operonEvents.setOperon2NegativeOrientation(operon2NegativeOrientation)
    operonEvents.setOperon1Reversed(operon1Reversed)
    operonEvents.setOperon2Reversed(operon2Reversed)
    
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