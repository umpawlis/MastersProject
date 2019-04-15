from SequenceService import computeOperonDifferences
from SequenceService import formatAllOperons
from SequenceService import findUniqueGenes
from SequenceService import addDuplicationEventsToStrain
from SequenceService import addDeletionEventsToStrain
from GenomeFragment import GenomeFragment
from Event import Event
import numpy as np
import globals
import copy

################################
###Global Alignment Functions###
################################

######################################################
# detectOrthologsByGlobalAlignment
# Parameters:
# Description: Calls function to construct the matrix and then calls function to scan matrix
######################################################
def findOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2):
    events = []
    globalAlignmentCounter = 0

    #Step 1: Compute the Matrices needed
    alignmentMatrix, eventMatrix = computeGlobalAlignmentMatrix(strain1, strain2)

    #Step 2: Scan and find all of the orthologous operons in the matrix
    events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = scanGlobalAlignmentMatrixForOrthologs(alignmentMatrix, eventMatrix, coverageTracker1, coverageTracker2, strain1, strain2)

    return events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2

######################################################
# computeGlobalAlignmentMatrix
# Parameters: directoryList - List of directories to process
# Description: Creates two matrices. a score matrix and a event matrix
######################################################
def computeGlobalAlignmentMatrix(strain1, strain2):

    firstOperonList = strain1.genomeFragments
    secondOperonList = strain2.genomeFragments
    strain1Name = strain1.name
    strain2Name = strain2.name

    print('Computing global alignment matrix for: {%s, %s}...' % (strain1Name, strain2Name))

    #initialize the matrix to store the global alignment scores
    globalAlignmentMatrix = [[ 0.0 for x in range(0, len(secondOperonList))] for y in range(0, len(firstOperonList))]
    eventMatrix = [[None for x in range(0, len(secondOperonList))] for y in range(0, len(firstOperonList))]

    ####################################
    ##Calculations
    ####################################
    for x in range(0, len(firstOperonList)):
        for y in range(0, len(secondOperonList)):
            op1 = firstOperonList[x] #Fragment
            op2 = secondOperonList[y] #Fragment

            #Case 1: An origin or terminus
            if op1.description == 'Origin' and op2.description == 'Origin' or op1.description == 'Terminus' and op2.description == 'Terminus':
                if op1.description == 'Origin':
                    print('Found Origin!')
                else:
                    print('Found Terminus!')

                event = Event(0)
                event.setScore(0)
                event.setGenome1Operon(op1.sequence)
                event.setGenome2Operon(op2.sequence)
                event.setGenome1Name(strain1Name)
                event.setGenome2Name(strain2Name)
                event.isOriginallyNegativeOrientationOp1(False)
                event.isOriginallyNegativeOrientationOp2(False)
                event.setOperon1Index(x)
                event.setOperon2Index(y)
                event.setTechnique(op1.description)
                event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                event.setOperon2Alignment(copy.deepcopy(op2.sequence))
                event.setFragmentDetails1(op1)
                event.setFragmentDetails2(op2)
                eventMatrix[x][y] = event
                globalAlignmentMatrix[x][y] = str(0) + '*'

            #Case 2: Two singleton genes are being compared
            elif len(op1.sequence) == 1 and len(op2.sequence) == 1:

                if op1.sequence[0] == op2.sequence[0]:
                    event = Event(0)
                    event.setScore(0)
                    event.setGenome1Operon(op1.sequence)
                    event.setGenome2Operon(op2.sequence)
                    event.setGenome1Name(strain1Name)
                    event.setGenome2Name(strain2Name)
                    event.isOriginallyNegativeOrientationOp1(op1.isNegativeOrientation)
                    event.isOriginallyNegativeOrientationOp2(op2.isNegativeOrientation)
                    event.setOperon1Index(x)
                    event.setOperon2Index(y)
                    event.setTechnique('Global Alignment')
                    event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                    event.setOperon2Alignment(copy.deepcopy(op2.sequence))
                    event.setFragmentDetails1(op1)
                    event.setFragmentDetails2(op2)
                    eventMatrix[x][y] = event
                    globalAlignmentMatrix[x][y] = str(0) + '*' #Singletons are a perfect match
                else:
                    globalAlignmentMatrix[x][y] = -999 #Singletons don't match

            #Case 3: Both are operons
            elif len(op1.sequence) > 1 and len(op2.sequence) > 1:
                event = Event(0)
                event.setGenome1Operon(op1.sequence)
                event.setGenome2Operon(op2.sequence)
                event.setGenome1Name(strain1Name)
                event.setGenome2Name(strain2Name)
                event.isOriginallyNegativeOrientationOp1(op1.isNegativeOrientation)
                event.isOriginallyNegativeOrientationOp2(op2.isNegativeOrientation)
                event.setOperon1Index(x)
                event.setOperon2Index(y)
                event.setTechnique('Global Alignment')
                event.setFragmentDetails1(op1)
                event.setFragmentDetails2(op2)

                score, event = performGlobalAlignment(op1.sequence, op2.sequence, event)

                event.setScore(score)
                eventMatrix[x][y] = event

                globalAlignmentMatrix[x][y] = score

                threshold = max(len(op1.sequence), len(op2.sequence))
                threshold = threshold//3
                numOperonDifferences = computeOperonDifferences(op1.sequence, op2.sequence)
                if numOperonDifferences <= threshold:
                    globalAlignmentMatrix[x][y] = str(globalAlignmentMatrix[x][y]) + '*'

            #Case 4: One of them is an operon and the other is a singleton
            else:
                globalAlignmentMatrix[x][y] = -999

    ####################################
    ##End of Calculations
    ####################################
    print ('Done computing global alignment matrix for {%s, %s}\n' % (strain1Name, strain2Name))
    #outputResultsToExcel(strain1Name, strain2Name, firstOperonList, secondOperonList, globalAlignmentMatrix)

    return globalAlignmentMatrix, eventMatrix

######################################################
# performGlobalAlignment
# Parameters:
# Description: Performs a global alignment
######################################################
def performGlobalAlignment(operon1, operon2, event):

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
                    scoreMatrix[a][b] = min(scoreMatrix[a-1][b-1] + globals.codonCost, scoreMatrix[a-1][b] + globals.deletionCost, scoreMatrix[a][b-1] + globals.deletionCost, scoreMatrix[a-1][b-1] + globals.substitutionCost)
            else:
                scoreMatrix[a][b] = min(scoreMatrix[a-1][b] + globals.deletionCost, scoreMatrix[a][b-1] + globals.deletionCost, scoreMatrix[a-1][b-1] + globals.substitutionCost)

    #Compute the number of events that occured between the operons
    event = globalAlignmentTraceback(scoreMatrix, operon1, operon2, event)

    return scoreMatrix[len(operon1)][len(operon2)], event

######################################################
# globalAlignmentTraceback
# Parameters:
# Description: Performs a traceback on a given matrix
######################################################
def globalAlignmentTraceback(matrix, operon1, operon2, event):

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

    #Track the alignment (two in the event we have substitutions)
    alignmentSequence1 = []
    alignmentSequence2 = []

    #Tracks where the extra genes are from
    gap1Indexes = []
    gap2Indexes = []

    while i > 0 or j > 0:
        #Case 1: Perfect match
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] and operon1[i-1] == operon2[j-1]:
            match += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + globals.codonCost) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            codonMismatch += 1

            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 3: Substitution
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + globals.substitutionCost):
            substitution += 1

            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and matrix[i][j] == (matrix[i-1][j] + globals.deletionCost):
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

    #Need to swap the gap lists since the gaps refer to extra genes
    temp = operon1Gaps
    operon1Gaps = operon2Gaps
    operon2Gaps = temp

    temp = gap1Indexes
    gap1Indexes = gap2Indexes
    gap2Indexes = temp

    event.setNumMatches(match)
    event.setNumGeneMismatches(mismatch)

    event.setNumCodonMismatches(codonMismatch)
    event.setNumSubstitutions(substitution)

    event.setOperon1Alignment(alignmentSequence1)
    event.setOperon2Alignment(alignmentSequence2)

    event.setOperon1Gaps(operon1Gaps)
    event.setOperon2Gaps(operon2Gaps)

    event.setOperon1GapIndexes(gap1Indexes)
    event.setOperon2GapIndexes(gap2Indexes)

    #Used for debugging
    #print('These are the operons being compared: %s, %s' %(operon1, operon2))
    #print('This is the resulting alignment: %s, %s' %(alignmentSequence1, alignmentSequence2))
    #print('These are the extra genes for operon 1: %s' %(operon1Gaps))
    #print('These are the indexes for extra genes in operon 1: %s' %(gap1Indexes))
    #print('These are the extra genes for operon 2: %s' %(operon2Gaps))
    #print('These are the indexes for extra genes in operon 2: %s' %(gap2Indexes))

    return event

######################################################
# scanGlobalAlignmentMatrixForOrthologs
# Parameters:
# Description: Scans matrix and identifies orthologous operons
######################################################
def scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, eventMatrix, coverageTracker1, coverageTracker2, strain1, strain2):
    events = []
    currentScoreSelected = 0
    globalAlignmentCounter = 0
    
    maxValue = findMaxValueInMatrix(globalAlignmentMatrix)

    #Keep iterating util we find all the optimal scores (Finding orthologs using global alignment)
    while currentScoreSelected <= maxValue:
        #Prioritize the selection of operons with the same sign
        for i in range(0, len(globalAlignmentMatrix)):
            for j in range(0, len(globalAlignmentMatrix[i])):
                #Check if this is a * score, if both operons have not been marked off and if both are the same orientation
                if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False) and (eventMatrix[i][j].originallyNegativeOrientationOp1 == eventMatrix[i][j].originallyNegativeOrientationOp2):
                    score = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                    #Check if the score matches the scores we're currently looking for
                    if score == currentScoreSelected:
                        #We found an ortholog in the global alignment matrix
                        print('\n##### Global Alignment #####')

                        globals.trackingId += 1
                        globalAlignmentCounter+=1

                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        event = eventMatrix[i][j]
                        event.trackingEventId = globals.trackingId
                        event, strain1, strain2 = reconstructOperonSequence(event, strain1, strain2)
                        event.printEvent()

                        #Add the event to the tracking events list
                        events.append(event)
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
                        globals.trackingId += 1
                        globalAlignmentCounter+=1

                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        event = eventMatrix[i][j]
                        event.trackingEventId = globals.trackingId
                        event, strain1, strain2 = reconstructOperonSequence(event, strain1, strain2)
                        event.printEvent()

                        #Add the event to the tracking events list
                        events.append(event)
                        print('###################################\n')
        currentScoreSelected += globals.codonCost

    return events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2

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
# reconstructOperonSequence
# Parameters:
# Description: Reconstructs the ancestral operon by determining whether the gaps are losses or duplications
######################################################
def reconstructOperonSequence(event, strain1, strain2):
    ancestralOperon = copy.deepcopy(event.operon1Alignment)
    if event.score == 0:
        print('No differences detected between these two operons')
        event.setAncestralOperonGeneSequence(ancestralOperon)
    else:
        print('Differences detected between these two operons!')
        operon1Gaps = event.operon1Gaps
        operon2Gaps = event.operon2Gaps
        operon1GapIndexes = event.operon1GapIndexes
        operon2GapIndexes = event.operon2GapIndexes

        print('These are the extra genes for operon 1: %s' %(operon1Gaps))
        print('These are the indexes for extra genes in operon 1: %s' %(operon1GapIndexes))
        print('These are the extra genes for operon 2: %s' %(operon2Gaps))
        print('These are the indexes for extra genes in operon 2: %s' %(operon2GapIndexes))

        #Checks if these extra genes are duplicates by checking if they exist within the alignment and removes them if they do
        operon1Gaps, duplicateSizesWithinAlignment1 = checkForMatchesInAlignment(operon1Gaps, event.operon1Alignment)
        operon2Gaps, duplicateSizesWithinAlignment2 = checkForMatchesInAlignment(operon2Gaps, event.operon2Alignment)

        #increment the duplicate counters
        #incrementDuplicateSizeCounters(duplicateSizesWithinAlignment1)
        #incrementDuplicateSizeCounters(duplicateSizesWithinAlignment2)
        strain1 = addDuplicationEventsToStrain(duplicateSizesWithinAlignment1, strain1)
        strain2 = addDuplicationEventsToStrain(duplicateSizesWithinAlignment2, strain2)

        i = len(operon1Gaps)
        j = len(operon2Gaps)

        while (i > 0) or (j > 0):
            #Select the gap with the biggest index b/c we will be performing the insertion rear to front of operon to avoid messing up the indexes of the other gaps
            if i > 0 and j > 0 and operon1GapIndexes[i-1] > operon2GapIndexes[j-1]:
                #This means both queues have gaps however the index in queue 1 is bigger so we'll deal with that one first
                #print('Gap being processed: %s' % (operon1Gaps[i]))
                numUniqueFound, deletionSizes, duplicationSizes, updateUnaligned = findUniqueGenes(operon1Gaps[i-1], strain1.formattedSequence, strain1.sequenceConversion[event.operon1Index])
                strain1 = addDuplicationEventsToStrain(duplicationSizes, strain1)
                strain2 = addDeletionEventsToStrain(deletionSizes, strain2)

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
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
                numUniqueFound, deletionSizes, duplicationSizes, updateUnaligned = findUniqueGenes(operon2Gaps[j-1], strain2.formattedSequence, strain2.sequenceConversion[event.operon2Index])
                strain2 = addDuplicationEventsToStrain(duplicationSizes, strain2)
                strain1 = addDeletionEventsToStrain(deletionSizes, strain1)

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
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
                numUniqueFound, deletionSizes, duplicationSizes, updateUnaligned = findUniqueGenes(operon1Gaps[i-1], strain1.formattedSequence, strain1.sequenceConversion[event.operon1Index])
                strain1 = addDuplicationEventsToStrain(duplicationSizes, strain1)
                strain2 = addDeletionEventsToStrain(deletionSizes, strain2)

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
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
                numUniqueFound, deletionSizes, duplicationSizes, updateUnaligned = findUniqueGenes(operon2Gaps[j-1], strain2.formattedSequence, strain2.sequenceConversion[event.operon2Index])
                strain2 = addDuplicationEventsToStrain(duplicationSizes, strain2)
                strain1 = addDeletionEventsToStrain(deletionSizes, strain1)

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
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
        event.setAncestralOperonGeneSequence(ancestralOperon)
        #print('This is the resulting ancestral operon: %s' % (ancestralOperon))
        #print('\n\n')
        #print('These are the extra genes remaining for operon 1: %s' %(operon1Gaps))
        #print('These are the extra genes remaining for operon 2: %s' %(operon2Gaps))
        #print('These are the duplicate sizes operon 1: %s' %(duplicateSizesWithinAlignment1))
        #print('These are the duplicate sizes operon 2: %s\n\n' %(duplicateSizesWithinAlignment2))

    return event, strain1, strain2

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