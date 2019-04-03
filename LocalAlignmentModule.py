from SequenceService import formatAllOperons
from SequenceService import reverseSequence
from SequenceService import formatAndComputeOperonDifferences
from SequenceService import findUniqueGenes
from SequenceService import incrementDuplicateSizeCounters
from SequenceService import incrementDeletionSizeCounters
from GlobalAlignmentModule import checkForMatchesInAlignment
from Event import Event
import globals
import numpy as np
import copy

################################
###Local Alignment Functions####
################################

######################################################
# detectOrthologsByLocalAlignment
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def findOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2):
    events = []
    localAlignmentCounter = 0
    minValue = min(len(coverageTracker1), len(coverageTracker2))
    sequence1 = strain1.getSequence()
    print "\nSequence:--\n"
    print sequence1
    sequence2 = strain2.getSequence()
    print sequence2
    genesStrain1 = strain1.getGenes()
    # print "\nGene Strain:-- \n"
    # print genesStrain1
    genesStrain2 = strain2.getGenes()
    operonPositionList1 = strain1.getOperonPositions()
    operonPositionList2 = strain2.getOperonPositions()
    singletonDict1 = strain1.getSingletonDict()
    singletonDict2 = strain2.getSingletonDict()
    formattedSequence1, operonIndexes1 = formatAllOperons(strain1.getSequence())
    # print "\nFormatted Sequence:--\n"
    # print formattedSequence1
    # print operonIndexes1
    formattedSequence2, operonIndexes2 = formatAllOperons(strain2.getSequence())
    
    #Scan the matrix x times to find optimal local alignments everytime an entire matrix is computed
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
            #Make sure the operon is not marked and not a singleton
            if coverageTracker1[i] == False and len(sequence1[i].split(',')) > 1:
                for j in range(0, len(coverageTracker2)):
                    #Make sure the second operon is not marked and not a singleton
                    if coverageTracker2[j] == False and len(sequence2[j].split(',')) > 1:
                        op1 = sequence1[i]
                        op2 = sequence2[j]
                        score, formattedOperon1, formattedOperon2, endPosition, startPosition, aligned1, aligned2, operonEvent = localAlignment(op1, op2, i, j, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)
                        
                        if (score > highestScore) or (score == highestScore and (abs(i - j)) < distance):
                            highestScore = score
                            rowIndex = i
                            colIndex = j
                            distance = abs(i - j)

                            # Testing delete next two line later
                            # chosenOp1 = op1
                            # chosenOp2 = op2

                            chosenOperon1 = formattedOperon1
                            chosenOperon2 = formattedOperon2
                            chosenStart = startPosition
                            chosenEnd = endPosition
                            chosenAligned1 = aligned1
                            chosenAligned2 = aligned2
                            chosenOpEvent = operonEvent

                            chosenOpEvent.setOperon1Index(i)
                            chosenOpEvent.setOperon2Index(j)

                            # printAlignments(i, j, formattedOperon1, formattedOperon2, aligned1, aligned2, "after extension")
        #After scanning the whole matrix if we found a best score, then store it
        if highestScore > -1:
            print('\n******** Local Alignment ***********')
            print('**************************************\n')

            # print "\nLa Operon 1: \n"
            # print chosenOp1
            # print rowIndex
            print operonPositionList1
            # print singletonDict1

            # print chosenOp2
            # print colIndex
            print operonPositionList2

            # print chosenStart
            # print chosenEnd

            globals.trackingId += 1
            localAlignmentCounter+=1
            coverageTracker1[rowIndex] = True
            coverageTracker2[colIndex] = True

            chosenOpEvent.setGenome1Name(strain1.getName())
            chosenOpEvent.setGenome2Name(strain2.getName())
            chosenOpEvent.setScore(highestScore)
            chosenOpEvent.trackingEventId = globals.trackingId

            ancestralOperon = determineAncestor(chosenOperon1, chosenOperon2, chosenStart, chosenEnd, chosenAligned1, chosenAligned2, formattedSequence1, formattedSequence2, operonIndexes1[rowIndex], operonIndexes2[colIndex], chosenOpEvent)
            chosenOpEvent.setAncestralOperonGeneSequence(ancestralOperon)
            chosenOpEvent.printEvent()
            # print "\n Operon Alignments \n"
            # print chosenOpEvent.operon1Alignment
            # print chosenOpEvent.operon2Alignment
            # print "\n Operon Gaps \n"
            # print chosenOpEvent.operon1Gaps
            # print chosenOpEvent.operon2Gaps
            # print "\n Operon Indexes \n"
            # print chosenOpEvent.operon1GapIndexes
            # print chosenOpEvent.operon2GapIndexes

            events.append(chosenOpEvent)

            #TODO: Create Event


            print('\n**************************************')
            print('**************************************\n\n')
                        
    return events, coverageTracker1, coverageTracker2, localAlignmentCounter

######################################################
# localAlignment
# Parameters:
# Description: Performs local alignment on two operons
######################################################
def localAlignment(op1, op2, op1Position, op2Position, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    #Local alignment parameters
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5
    mismatch = -1.0
    gap = -1.0
    
    leftAdjustment1 = 0
    leftAdjustment2 = 0
    rightAdjustment1 = 0
    rightAdjustment2 = 0
    
    #Check whether the operons are in the - orientation
    negativeOrientationOp1 = reverseSequence(op1)
    negativeOrientationOp2 = reverseSequence(op2)
    
    #Tracks whether we reversed the operons
    operon1Reversed = False
    operon2Reversed = False
    
    #Format and compute the operon differences
    setDifference, operon1, operon2, numberOfDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

    #Create event for this comparison
    event = Event(0)
    event.setGenome1Operon(operon1)
    event.setGenome2Operon(operon2)
    event.isOriginallyNegativeOrientationOp1(negativeOrientationOp1) #This tracks the original orientation of op1
    event.isOriginallyNegativeOrientationOp2(negativeOrientationOp2) #This tracks the original orientation of op2
    event.setTechnique('Local Alignment')
    
    #If the orientation of the operons does not match, then flip the operon in the negative orientation to the positive orientation
    if negativeOrientationOp1 != negativeOrientationOp2:
        if negativeOrientationOp1:
            operon1.reverse()
            operon1Reversed = True
            negativeOrientationOp1 = False
            
        if negativeOrientationOp2:
            operon2.reverse()
            operon2Reversed = True
            negativeOrientationOp2 = False
            
    if negativeOrientationOp1 != negativeOrientationOp2:
        print('Check code! These operons should be in the same orientation!')

    #Track whether these operons were reversed
    event.isReversedOp1(operon1Reversed) #This tracks whether op1 was reversed
    event.isReversedOp2(operon2Reversed) #This tracks whether op2 was reversed
        
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
                #Changed the code a little bit to make sure the opmital choice is made (codon mismatches were not always optimal in global alignments)
                else:
                    scoreMatrix[a][b] = max(scoreMatrix[a-1][b-1] + matchWithoutCodon, scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)
            else:
                scoreMatrix[a][b] = max(scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)

            if scoreMatrix[a][b] > maxScore:
                maxScore = scoreMatrix[a][b]
                maxPosition = (a, b)
                
    #Trace back score matrix to get alignment
    aligned1, aligned2, numGaps, endPosition = traceback(operon1, operon2, scoreMatrix, maxPosition, event)
    # print aligned1
    # print aligned2

    shortestLength = min(len(operon1), len(operon2))
    
    #Adjust gene position we are looking at for extension based on the orienation of the operon
    if negativeOrientationOp1:
        leftAdjustment1 = len(operon1)-1
        leftAdjustment1 -= endPosition[0]-1
        rightAdjustment1 += (len(operon1) - (maxPosition[0]))
    else:
        leftAdjustment1 += endPosition[0]-1
        rightAdjustment1 = len(operon1)-1
        rightAdjustment1 -= (len(operon1) - (maxPosition[0]))

    if negativeOrientationOp2:
        leftAdjustment2 = len(operon2)-1
        leftAdjustment2 -= endPosition[1]-1
        rightAdjustment2 += (len(operon2) - (maxPosition[1]))
    else:
        leftAdjustment2 += endPosition[1]-1
        rightAdjustment2 = len(operon2)-1
        rightAdjustment2 -= (len(operon2) - (maxPosition[1]))

    #Only extend operons when the current alignment has no gaps
    # if genesStrain1 and genesStrain2:
    #     if len(aligned1) == shortestLength:
    #         #print("One operon is a SUBSET of the other. Trying extention..")
    #         extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, negativeOrientationOp1, negativeOrientationOp2, aligned1, aligned2, singletonDict1, singletonDict2)
    #         returningScore = maxScore + extensionScore
    #         extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, negativeOrientationOp1, negativeOrientationOp2, aligned1, aligned2, singletonDict1, singletonDict2)
    #         returningScore += extensionScore
    #         printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left and right extension")
    #     elif len(aligned1) < shortestLength:
    #         if maxPosition[0] == len(scoreMatrix)-1 or maxPosition[1] == len(scoreMatrix[0])-1:
    #             #print("Trying right extension..")
    #             extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, negativeOrientationOp1, negativeOrientationOp2, aligned1, aligned2, singletonDict1, singletonDict2)
    #             if len(aligned1) >= shortestLength:
    #                 printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after right extension")
    #                 returningScore = maxScore + extensionScore
    #         elif endPosition[0] == 1 or endPosition[1] == 1:
    #             #print("Trying left extension..")
    #             extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, negativeOrientationOp1, negativeOrientationOp2, aligned1, aligned2, singletonDict1, singletonDict2)
    #             if len(aligned1) >= shortestLength:
    #                 printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left extension")
    #                 returningScore = maxScore + extensionScore
    #No extension but still qualifies because it has no gaps in alignment
    if len(aligned1) >= (0.75 * shortestLength):
        printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after no extension due to missing gene list(s)")
        returningScore = maxScore
    #print("Final alignment:")
    #print(aligned1)
    #print(aligned2)

    return returningScore, operon1, operon2, maxPosition, endPosition, aligned1, aligned2, event

######################################################
# traceback
# Parameters: operon1, operon2, scoreMatrix, startPosition
# Description: treaverses score matrix to determine optimal alignment
######################################################
def traceback(operon1, operon2, scoreMatrix, startPosition, event):
    match = 1
    codonMismatch = 0
    mismatch = 0

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = startPosition
    move = nextMove(scoreMatrix, x, y)
    numGaps = 0

    op1Gaps = []
    op1Gap = []
    op1ConsecutiveGap = False #Tracks consecutive gaps

    op2Gaps = []
    op2Gap = []
    op2ConsecutiveGap = False #Tracks consecutive gaps

    #Tracks where the extra genes are from
    gap1Indexes = []
    gap2Indexes = []

    while move != END:
        if move == DIAG:
            if operon1[x-1] == operon2[y-1]:
                match += 1
            elif operon1[x-1].split('_')[0].strip() == operon2[y-1].split('_')[0].strip():
                codonMismatch += 1
            else:
                mismatch += 1

            aligned_seq1.insert(0, operon1[x - 1])
            aligned_seq2.insert(0, operon2[y - 1])
            x -= 1
            y -= 1

        elif move == UP:
            index = x - 1
            # aligned_seq1.insert(0, operon1[x - 1])
            # aligned_seq2.insert(0, ' - ')
            x -= 1
            numGaps += 1

            op1ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if op2ConsecutiveGap:
                op2Gap.insert(0, operon1[index])
                op2ConsecutiveGap = True
            else:
                if len(op2Gap) > 0:
                    op2Gaps.insert(0, op2Gap)
                op2Gap = []
                op2Gap.insert(0, operon1[index])
                gap2Indexes.insert(0, len(aligned_seq2))
                op2ConsecutiveGap = True

        else:
            index = y - 1
            # aligned_seq1.insert(0, ' - ')
            # aligned_seq2.insert(0, operon2[y - 1])
            y -= 1
            numGaps += 1

            op2ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if op1ConsecutiveGap:
                op1Gap.insert(0, operon2[index])
                op1ConsecutiveGap = True
            else:
                if len(op1Gap) > 0:
                    op1Gaps.insert(0, op1Gap)
                op1Gap = []
                op1Gap.insert(0, operon2[index])
                gap1Indexes.insert(0, len(aligned_seq1))
                op1ConsecutiveGap = True

        move = nextMove(scoreMatrix, x, y)

    aligned_seq1.insert(0, operon1[x - 1])
    aligned_seq2.insert(0, operon2[y - 1])
    endPosition = (x,y)

    #Empty any remaining gaps
    if len(op1Gap) > 0:
        op1Gaps.insert(0, op1Gap)
    if len(op2Gap) > 0:
        op2Gaps.insert(0, op2Gap)

    #The indexes values need to be flipped b/c right now they're oriented from right to left
    if len(gap1Indexes) > 0:
        for x in range(0, len(gap1Indexes)):
            gap1Indexes[x] = len(aligned_seq1) - gap1Indexes[x]
    if len(gap2Indexes) > 0:
        for x in range(0, len(gap2Indexes)):
            gap2Indexes[x] = len(aligned_seq2) - gap2Indexes[x]

    event.setNumMatches(match)
    event.setNumCodonMismatches(codonMismatch)
    event.setNumGeneMismatches(mismatch)
    # event.setNumSubstitutions(substitution)
    # print " Operons: "
    # print operon1
    # print operon2
    # print " Alignment: "
    # print aligned_seq1
    # print aligned_seq2
    event.setOperon1Alignment(aligned_seq1)
    event.setOperon2Alignment(aligned_seq2)
    # The gaps we found in op2, we will use to look for duplicates/losses in op1, and vice versa
    event.setOperon1Gaps(op2Gaps)
    event.setOperon2Gaps(op1Gaps)
    event.setOperon1GapIndexes(gap2Indexes)
    event.setOperon2GapIndexes(gap1Indexes)

    return aligned_seq1, aligned_seq2, numGaps, endPosition

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

    # print ("Performing " + direction + " extension: ")
    # print operon1
    # print operon2
    # print ("With following alignment ")
    # print aligned1
    # print aligned2
    # print ("Gene Positions")
    # print opGenePosition1
    # print opGenePosition2

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
                    # print("Both singletons")
            elif (str(opGenePosition1) in singletonDict1 and opGenePosition2 in operonRange2):
                if (('-' in singletonDict1[str(opGenePosition1)]) != reverseOp1):
                    mismatch = True
                    # print("gene1 is singletons")
            elif (str(opGenePosition2) in singletonDict2 and opGenePosition1 in operonRange1):
                if (('-' in singletonDict2[str(opGenePosition2)]) != reverseOp2):
                    mismatch = True
                    # print("gene2 is singletons")
            else:
                mismatch = True
                # print ("One or both in a separate operon")

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
                    # print ("Inserting genes:")
                    # print aligned1
                    # print aligned2
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
def determineAncestor(op1, op2, startPosition, endPosition, aligned1, aligned2, sequence1, sequence2, opIndex1, opIndex2, event):
    ancestralOperon = copy.deepcopy(event.operon1Alignment)
    updatedUnaligned = []
    unaligned1 = []
    unaligned2 = []
    deletionSizes = None
    duplicationSizes = None

    print "Aligned Sections"
    print aligned1
    print aligned2


    unaligned = getUnaligned(op1, startPosition[0]-1, endPosition[0]-1)
    print "Unaligned sections"
    print unaligned
    #numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned2, unaligned)
    numUnique = len(op1)-len(aligned2)
    numDuplicateFound = 0
    alignedRange = (startPosition[0]-1, endPosition[0]-1)

    for geneList in unaligned:
        if geneList:
            numUniqueFound, deletionSizes, duplicationSizes, updatedGeneList = findUniqueGenes(geneList, sequence1, opIndex1, alignedRange)
            numUnique -= numUniqueFound
            numDuplicateFound += numUniqueFound
            updatedUnaligned.append(updatedGeneList)
        else:
            updatedUnaligned.append([])

    unaligned = getUnaligned(op2, startPosition[1]-1, endPosition[1]-1)
    print "Unaligned sections"
    print unaligned
    # numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned1, unaligned)
    numUnique = len(op2)-len(aligned1)
    numDuplicateFound = 0
    alignedRange = (startPosition[1]-1, endPosition[1]-1)

    for geneList in unaligned:
        if geneList:
            numUniqueFound, deletionSizes, duplicationSizes, updatedGeneList = findUniqueGenes(geneList, sequence2, opIndex2, alignedRange)
            numUnique -= numUniqueFound
            numDuplicateFound += numUniqueFound
            updatedUnaligned.append(updatedGeneList)
        else:
            updatedUnaligned.append([])

    print updatedUnaligned

    while updatedUnaligned[0] or updatedUnaligned[2]:
        if updatedUnaligned[0]:
            gene = updatedUnaligned[0].pop()
            if gene != '-':
                unaligned1.insert(0, gene)
        if updatedUnaligned[2]:
            gene = updatedUnaligned[2].pop()
            if gene != '-':
                unaligned1.insert(0, gene)

    while updatedUnaligned[1] or updatedUnaligned[3]:
        if updatedUnaligned[1]:
            gene = updatedUnaligned[1].pop()
            if gene != '-':
                unaligned2.insert(0, gene)
        if updatedUnaligned[3]:
            gene = updatedUnaligned[3].pop()
            if gene != '-':
                unaligned2.insert(0, gene)

    print('Differences detected between these two operons!')
    operon1Gaps = event.operon1Gaps
    operon2Gaps = event.operon2Gaps
    operon1GapIndexes = event.operon1GapIndexes
    operon2GapIndexes = event.operon2GapIndexes

    print operon1Gaps
    print operon1GapIndexes
    print operon2Gaps
    print operon2GapIndexes

    print('These are the extra genes for operon 1: %s' %(operon1Gaps))
    print('These are the indexes for extra genes in operon 1: %s' %(operon1GapIndexes))
    print('These are the extra genes for operon 2: %s' %(operon2Gaps))
    print('These are the indexes for extra genes in operon 2: %s' %(operon2GapIndexes))

    #Checks if these extra genes are duplicates by checking if they exist within the alignment and removes them if they do
    operon1Gaps, duplicateSizesWithinAlignment1 = checkForMatchesInAlignment(operon1Gaps, event.operon1Alignment)
    operon2Gaps, duplicateSizesWithinAlignment2 = checkForMatchesInAlignment(operon2Gaps, event.operon2Alignment)
    
    #increment the duplicate counters
    incrementDuplicateSizeCounters(duplicateSizesWithinAlignment1)
    incrementDuplicateSizeCounters(duplicateSizesWithinAlignment2)
    
    i = len(operon1Gaps)
    j = len(operon2Gaps)

    while (i > 0) or (j > 0):
        #Select the gap with the biggest index b/c we will be performing the insertion rear to front of operon to avoid messing up the indexes of the other gaps
        if i > 0 and j > 0 and operon1GapIndexes[i-1] > operon2GapIndexes[j-1]:
            #This means both queues have gaps however the index in queue 1 is bigger so we'll deal with that one first
            #print('Gap being processed: %s' % (operon1Gaps[i]))
            numUniqueFound, deletionSizes, duplicationSizes, _ = findUniqueGenes(operon1Gaps[i-1], sequence1, opIndex1)
            incrementDuplicateSizeCounters(duplicationSizes)
            incrementDeletionSizeCounters(deletionSizes)
            #print('Gap being processed: %s' % (operon1Gaps[i-1]))
            #print('Number of unique genes found: %s' %(numUniqueFound))
            #print('Number of deletion genes found: %s' %(deletionSizes))
            #print('Number of duplicate genes found: %s' %(duplicationSizes))
            if len(operon1Gaps[i-1]) > 0:
                #Insert gap into operon
                operon1Gaps[i-1].reverse()
                for gene in operon1Gaps[i-1]:
                    ancestralOperon.insert(operon1GapIndexes[i-1], gene)
            i = i - 1
        elif i > 0 and j > 0 and operon1GapIndexes[i-1] < operon2GapIndexes[j-1]:
            #This means both queues have gaps however the index in queue 2 is bigger so we'll insert that one first
            #print('Gap being processed: %s' % (operon2Gaps[j-1]))
            numUniqueFound, deletionSizes, duplicationSizes, _ = findUniqueGenes(operon2Gaps[j-1], sequence2, opIndex2)
            incrementDuplicateSizeCounters(duplicationSizes)
            incrementDeletionSizeCounters(deletionSizes)
            #print('Gap being processed: %s' % (operon2Gaps[j-1]))
            #print('Number of unique genes found: %s' %(numUniqueFound))
            #print('Number of deletion genes found: %s' %(deletionSizes))
            #print('Number of duplicate genes found: %s' %(duplicationSizes))
            if len(operon2Gaps[j-1]) > 0:
                #Insert gap into operon
                operon2Gaps[j-1].reverse()
                for gene in operon2Gaps[j-1]:
                    ancestralOperon.insert(operon2GapIndexes[j-1], gene)
            j = j - 1
        elif i > 0:
            #This means that queue 2 has no more gaps so we process the remaining gaps in queue 1
            #print('Gap being processed: %s' % (operon1Gaps[i-1]))
            numUniqueFound, deletionSizes, duplicationSizes, _ = findUniqueGenes(operon1Gaps[i-1], sequence1, opIndex1)
            incrementDuplicateSizeCounters(duplicationSizes)
            incrementDeletionSizeCounters(deletionSizes)
            #print('Gap being processed: %s' % (operon1Gaps[i-1]))
            #print('Number of unique genes found: %s' %(numUniqueFound))
            #print('Number of deletion genes found: %s' %(deletionSizes))
            #print('Number of duplicate genes found: %s' %(duplicationSizes))
            if len(operon1Gaps[i-1]) > 0:
                #Insert gap into operon
                operon1Gaps[i-1].reverse()
                for gene in operon1Gaps[i-1]:
                    ancestralOperon.insert(operon1GapIndexes[i-1], gene)
            i = i - 1
        elif j > 0:
            #This means that queue 1 has no more gaps to process so we deal with the remaining gaps in queue 2
            #print('Gap being processed: %s' % (operon2Gaps[j-1]))
            numUniqueFound, deletionSizes, duplicationSizes, _ = findUniqueGenes(operon2Gaps[j-1], sequence2, opIndex2)
            incrementDuplicateSizeCounters(duplicationSizes)
            incrementDeletionSizeCounters(deletionSizes)
            #print('Gap being processed: %s' % (operon2Gaps[j-1]))
            #print('Number of unique genes found: %s' %(numUniqueFound))
            #print('Number of deletion genes found: %s' %(deletionSizes))
            #print('Number of duplicate genes found: %s' %(duplicationSizes))
            if len(operon2Gaps[j-1]) > 0:
                #Insert gap into operon
                operon2Gaps[j-1].reverse()
                for gene in operon2Gaps[j-1]:
                    ancestralOperon.insert(operon2GapIndexes[j-1], gene)
            j = j - 1
    #Set ancestral operon
    ancestralOperon = unaligned1 + ancestralOperon + unaligned2
    event.setAncestralOperonGeneSequence(ancestralOperon)

    return ancestralOperon

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
# printStrains, printAlignments, printStats
# Parameters:
# Description: Prints a summary of the local alignments.
######################################################
def printAlignments(op1Position, op2Position, operon1, operon2, alignment1, alignment2, message):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("Strain 1 Operon %s: %s\n" %(op1Position, operon1))
    file.write("Strain 2 Operon  %s: %s\n" %(op2Position, operon2))
    file.write("Alignment Result (%s):\n" %message)
    file.write("%s\n" %alignment1)
    file.write("%s\n\n" %alignment2)

    file.close()