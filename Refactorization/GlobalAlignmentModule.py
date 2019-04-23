from SequenceService import addDuplicationEventsToStrain
from SequenceService import addDeletionEventsToStrain
from SequenceService import computeOperonDifferences
from Event import Event
import numpy as np
import globals
import copy

##################################
### Global Alignment Functions ###
##################################

######################################################
# findOrthologsByGlobalAlignment
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

    print('Computing global alignment matrix for: {%s, %s}...' % (strain1.name, strain2.name))

    #initialize the matrix to store the global alignment scores
    globalAlignmentMatrix = [[ 0.0 for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))]
    eventMatrix = [[None for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))]

    ####################################
    ##Calculations
    ####################################
    for x in range(0, len(strain1.genomeFragments)):
        for y in range(0, len(strain2.genomeFragments)):
            op1 = strain1.genomeFragments[x] #Fragment
            op2 = strain2.genomeFragments[y] #Fragment

            event = Event(0)
            event.setScore(0.0)
            event.setFragmentDetails1(op1)
            event.setFragmentDetails2(op2)
            event.setGenome1Name(strain1.name)
            event.setGenome2Name(strain2.name)
            event.setTechnique('Global Alignment')

            #Case 1: An origin or terminus
            if (op1.description == 'Origin' and op2.description == 'Origin') or (op1.description == 'Terminus' and op2.description == 'Terminus'):
                event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                event.setOperon2Alignment(copy.deepcopy(op2.sequence))

                eventMatrix[x][y] = event
                globalAlignmentMatrix[x][y] = str(0) + '*'

            #Case 2: Two singleton genes are being compared
            elif len(op1.sequence) == 1 and len(op2.sequence) == 1:

                if op1.sequence[0] == op2.sequence[0]:
                    event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                    event.setOperon2Alignment(copy.deepcopy(op2.sequence))

                    eventMatrix[x][y] = event
                    globalAlignmentMatrix[x][y] = str(0) + '*'
                else:
                    globalAlignmentMatrix[x][y] = -999 #Singletons don't match

            #Case 3: Both are operons
            elif len(op1.sequence) > 1 and len(op2.sequence) > 1:
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
    print ('Done computing global alignment matrix for {%s, %s}\n' % (strain1.name, strain2.name))
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

    #Tracks index and genes for codon mismatches in both strains
    codonMismatchIndexesStrain1 = []
    codonMismatchIndexesStrain2 = []
    codonMismatchGenesStrain1 = []
    codonMismatchGenesStrain2 = []

    #Tracks substitution indexes and genes for both strains
    substitutionIndexesStrain1 = []
    substitutionIndexesStrain2 = []
    substitutionGenesStrain1 = []
    substitutionGenesStrain2 = []

    #Tracks the genes in a gap and the index of those genes in both strains note: details stored in an array of arrays
    operon1Gaps = []
    operon1Gap = []
    operon1GapIndexes = [] #This is used to determine the position of the genes with respect to the genome
    operon1GapIndex = []
    operon1ConsecutiveGap = False #Tracks consecutive gaps

    operon2Gaps = []
    operon2Gap = []
    operon2GapIndexes = []
    operon2GapIndex = []
    operon2ConsecutiveGap = False #Tracks consecutive gaps

    #Track the alignment (two in the event we have substitutions)
    alignmentSequence1 = []
    alignmentSequence2 = []

    #Tracks where the extra genes are from
    gap1Indexes = [] #This is used to determine where to insert the genes into the alignment
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

            codonMismatchIndexesStrain1.append(i-1)
            codonMismatchGenesStrain1.append(operon1[i-1])

            codonMismatchIndexesStrain2.append(j-1)
            codonMismatchGenesStrain2.append(operon2[j-1])

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
        #Case 3: Substitution
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + globals.substitutionCost):
            substitution += 1

            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            substitutionIndexesStrain1.append(i-1)
            substitutionGenesStrain1.append(operon1[i-1])

            substitutionIndexesStrain2.append(j-1)
            substitutionGenesStrain2.append(operon2[j-1])

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
                operon2GapIndex.insert(0, index)

                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.insert(0, operon2Gap)
                    operon2GapIndexes.insert(0, operon2GapIndex)
                operon2Gap = []
                operon2GapIndex = []
                operon2Gap.insert(0, operon1[index])
                operon2GapIndex.insert(0, index)
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
                operon1GapIndex.insert(0, index)

                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.insert(0, operon1Gap)
                    operon1GapIndexes.insert(0, operon1GapIndex)
                operon1Gap = []
                operon1GapIndex = []
                operon1Gap.insert(0, operon2[index])
                operon1GapIndex.insert(0, index)
                gap1Indexes.insert(0, len(alignmentSequence1))
                operon1ConsecutiveGap = True

    #Empty any remaining gaps
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1GapIndexes.insert(0, operon1GapIndex)
        operon1Gap = []
        operon1GapIndex = []

    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2GapIndexes.insert(0, operon2GapIndex)
        operon2Gap = []
        operon2GapIndex = []

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

    temp = operon1GapIndexes
    operon1GapIndexes = operon2GapIndexes
    operon2GapIndexes = temp

    temp = gap1Indexes
    gap1Indexes = gap2Indexes
    gap2Indexes = temp

    #Match, mismatch details
    event.setNumMatches(match)
    event.setNumMismatches(mismatch)

    #Codon Mismatch details
    event.setNumCodonMismatches(codonMismatch)
    event.setCodonMismatchGenesStrain1(codonMismatchGenesStrain1)
    event.setCodonMismatchGenesStrain2(codonMismatchGenesStrain2)
    event.setCodonMismatchIndexesStrain1(codonMismatchIndexesStrain1)
    event.setCodonMismatchIndexesStrain2(codonMismatchIndexesStrain2)

    #Substitution details
    event.setNumSubstitutions(substitution)
    event.setSubstitutionIndexesStrain1(substitutionIndexesStrain1)
    event.setSubstitutionIndexesStrain2(substitutionIndexesStrain2)
    event.setSubstitutionGenesStrain1(substitutionGenesStrain1)
    event.setSubstitutionGenesStrain2(substitutionGenesStrain2)

    #Alignment details
    event.setOperon1Alignment(alignmentSequence1)
    event.setOperon2Alignment(alignmentSequence2)

    #Gap details
    event.setOperon1Gaps(operon1Gaps)
    event.setOperon2Gaps(operon2Gaps)
    event.setOperon1GapPositions(operon1GapIndexes)
    event.setOperon2GapPositions(operon2GapIndexes)
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
def scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, eventMatrix, coverageTracker1, coverageTracker2, strain1, strain2):
    events = []
    currentScoreSelected = 0
    globalAlignmentCounter = 0

    maxValue = findMaxValueInMatrix(globalAlignmentMatrix)

    #Keep iterating util we find all the optimal scores (Finding orthologs using global alignment)
    while currentScoreSelected <= maxValue:
        #For each cell in the matrix
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

                        #Codon mismatches
                        if event.numCodonMismatches > 0:
                            codonMismatchDescription1 = constructStatement(event.codonMismatchIndexesStrain1, event.codonMismatchGenesStrain1, event.fragmentDetails1)
                            codonMismatchDescription2 = constructStatement(event.codonMismatchIndexesStrain2, event.codonMismatchGenesStrain2, event.fragmentDetails2)

                            strain1.addCodonMismatchDetails(codonMismatchDescription1)
                            strain1.addCodonMismatchDetails(codonMismatchDescription2)

                        #Substitutions
                        if event.numSubstitutions > 0:
                            substitutionDescription1 = constructStatement(event.setSubstitutionIndexesStrain1, event.setSubstitutionGenesStrain1, event.fragmentDetails1)
                            substitutionDescription2 = constructStatement(event.setSubstitutionIndexesStrain2, event.setSubstitutionGenesStrain2, event.fragmentDetails2)

                            strain1.addSubstitutionDetails(substitutionDescription1)
                            strain1.addSubstitutionDetails(substitutionDescription2)

                        print(event.toString())

                        #Add the event to the tracking events list
                        events.append(event)
                        print('###################################\n')

        currentScoreSelected += globals.codonCost
    return events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2

######################################################
# constructStatement
# Parameters:
# Description: constructs the details of where the codon mismatch or substitution occured in the geneome
######################################################
def constructStatement(indexes, genes, fragmentDetails):
    temp = ''
    #TODO keep in mind the scenario if the operon was in the negative orientation originally
    startPosition = fragmentDetails.startPositionInGenome
    if len(indexes) != len(genes):
        print('Error! These two arrays should be the same length for Codon Mismatch Substitution parallel arrays')
    else:
        for x in range(0, len(indexes)):
            position = indexes[x] + startPosition
            gene = genes[x]
            temp+=gene + ' ' + str(position) + ';'

    return temp

######################################################
# reconstructOperonSequence
# Parameters:
# Description: Reconstructs the ancestral operon by determining whether the gaps are losses or duplications
######################################################
def reconstructOperonSequence(event, strain1, strain2):
    ancestralOperon = copy.deepcopy(event.operon1Alignment) #The alignments will be identical except when codon mismatches or substitutions occur
    if event.score == 0:
        print('No differences detected between these two operons')
        event.setAncestralOperonGeneSequence(ancestralOperon)
    else:
        print('Differences detected between these two operons!')
        operon1Gaps = event.operon1Gaps
        operon2Gaps = event.operon2Gaps
        operon1GapIndexes = event.operon1GapIndexes
        operon2GapIndexes = event.operon2GapIndexes
        operon1GapPositions = event.operon1GapPositions
        operon2GapPositions = event.operon2GapPositions

        print('These are the extra genes for operon 1: %s' %(operon1Gaps))
        print('These are the indexes for extra genes in operon 1: %s' %(operon1GapIndexes))
        print('These are the positions of the extra genes in operon 1: %s' %(operon1GapPositions))
        print('These are the extra genes for operon 2: %s' %(operon2Gaps))
        print('These are the indexes for extra genes in operon 2: %s' %(operon2GapIndexes))
        print('These are the positions of the extra genes in operon 2: %s' %(operon2GapPositions))

        #Step 1: Check if these extra genes are the result of a duplicate event within the alignment, remove them if they are
        operon1Gaps, operon1GapPositions, duplicateSizesWithinAlignment1, duplicationDetails1 = checkForMatchesWithinAlignment(operon1Gaps, event.operon1Alignment, operon1GapPositions, event.fragmentDetails1)
        operon2Gaps, operon2GapPositions, duplicateSizesWithinAlignment2, duplicationDetails2 = checkForMatchesWithinAlignment(operon2Gaps, event.operon2Alignment, operon2GapPositions, event.fragmentDetails2)

        #Add the details to the respective strain
        strain1 = addDuplicationEventsToStrain(strain1, duplicateSizesWithinAlignment1, duplicationDetails1)
        strain1 = addDuplicationEventsToStrain(strain1, duplicateSizesWithinAlignment1, duplicationDetails1)

        #Step 2: Check if these extra genes are the result of a duplication event within another operon, remove them if they are, else insert them
        i = len(operon1Gaps)
        j = len(operon2Gaps)

        while (i > 0) or (j > 0):
            #Select the gap with the biggest index b/c we will be performing the insertion rear to front of operon to avoid messing up the indexes of the other gaps
            if i > 0 and j > 0 and operon1GapIndexes[i-1] > operon2GapIndexes[j-1]:
                #This means both queues have gaps however the index in queue 1 is bigger so we'll deal with that one first
                #print('Gap being processed: %s' % (operon1Gaps[i]))
                duplicationSizes, duplicationDetails, operon1Gaps[i-1], operon1GapPositions[i-1] = checkForMatchesWithinOperons(strain1.genomeFragments, event.fragmentDetails1, operon1Gaps[i-1], operon1GapPositions[i-1])
                strain1 = addDuplicationEventsToStrain(strain1, duplicationSizes, duplicationDetails) #Adds duplication details to strain

                #print('Gap being processed: %s' % (operon1Gaps[i-1]))
                #print('Number of unique genes found: %s' %(numUniqueFound))
                #print('Number of deletion genes found: %s' %(deletionSizes))
                #print('Number of duplicate genes found: %s' %(duplicationSizes))
                if len(operon1Gaps[i-1]) > 0:
                    deletionDetails = ''
                    deletionSizes = []
                    #Insert gap into operon
                    operon1Gaps[i-1].reverse()
                    for k in range(0, len(operon1Gaps[i-1])):
                        ancestralOperon.insert(operon1GapIndexes[i-1], operon1Gaps[i-1][k])
                        deletionDetails += operon1Gaps[i-1][k] + ' ' + str(operon1GapPositions[i-1][k] + event.fragmentDetails2.startPositionInGenome) + ',' #TODO this might have to be fixed if the operon is in a negative orientation
                    deletionDetails += '|'                      #End of deleted segment
                    deletionSizes.append(len(operon1Gaps[i-1])) #Size of segment
                    strain1 = addDeletionEventsToStrain(strain1, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
                i = i - 1
            elif i > 0 and j > 0 and operon1GapIndexes[i-1] < operon2GapIndexes[j-1]:
                #This means both queues have gaps however the index in queue 2 is bigger so we'll insert that one first
                #print('Gap being processed: %s' % (operon2Gaps[j-1]))
                duplicationSizes, duplicationDetails, operon2Gaps[j-1], operon2GapPositions[j-1] = checkForMatchesWithinOperons(strain2.genomeFragments, event.fragmentDetails2, operon2Gaps[j-1], operon2GapPositions[j-1])
                strain2 = addDuplicationEventsToStrain(strain2, duplicationSizes, duplicationDetails) #Adds duplication details to strain

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
                #print('Gap being processed: %s' % (operon2Gaps[j-1]))
                #print('Number of unique genes found: %s' %(numUniqueFound))
                #print('Number of deletion genes found: %s' %(deletionSizes))
                #print('Number of duplicate genes found: %s' %(duplicationSizes))
                if len(operon2Gaps[j-1]) > 0:
                    deletionDetails = ''
                    deletionSizes = []
                    #Insert gap into operon
                    operon2Gaps[j-1].reverse()
                    for k in range (0, len(operon2Gaps[j-1])):
                        ancestralOperon.insert(operon2GapIndexes[j-1], operon2Gaps[j-1][k])
                        deletionDetails += operon2Gaps[j-1][k] + ' ' + str(operon2GapPositions[j-1][k] + event.fragmentDetails1.startPositionInGenome) + ',' #TODO this might change if operon is in negative orientation
                    deletionDetails += '|'                      #End of deleted segment
                    deletionSizes.append(len(operon2Gaps[j-1])) #Size of segment
                    strain2 = addDeletionEventsToStrain(strain2, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
                j = j - 1
            elif i > 0:
                #This means that queue 2 has no more gaps so we process the remaining gaps in queue 1
                #print('Gap being processed: %s' % (operon1Gaps[i-1]))
                duplicationSizes, duplicationDetails, operon1Gaps[i-1], operon1GapPositions[i-1] = checkForMatchesWithinOperons(strain1.genomeFragments, event.fragmentDetails1, operon1Gaps[i-1], operon1GapPositions[i-1])
                strain1 = addDuplicationEventsToStrain(strain1, duplicationSizes, duplicationDetails) #Adds duplication details to strain

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
                #print('Gap being processed: %s' % (operon1Gaps[i-1]))
                #print('Number of unique genes found: %s' %(numUniqueFound))
                #print('Number of deletion genes found: %s' %(deletionSizes))
                #print('Number of duplicate genes found: %s' %(duplicationSizes))
                if len(operon1Gaps[i-1]) > 0:
                    deletionDetails = ''
                    deletionSizes = []
                    #Insert gap into operon
                    operon1Gaps[i-1].reverse()
                    for k in range(0, len(operon1Gaps[i-1])):
                        ancestralOperon.insert(operon1GapIndexes[i-1], operon1Gaps[i-1][k])
                        deletionDetails += operon1Gaps[i-1][k] + ' ' + str(operon1GapPositions[i-1][k] + event.fragmentDetails2.startPositionInGenome) + ',' #TODO this might have to be fixed if the operon is in a negative orientation
                    deletionDetails += '|'                      #End of deleted segment
                    deletionSizes.append(len(operon1Gaps[i-1])) #Size of segment
                    strain1 = addDeletionEventsToStrain(strain1, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
                i = i - 1
            elif j > 0:
                #This means that queue 1 has no more gaps to process so we deal with the remaining gaps in queue 2
                #print('Gap being processed: %s' % (operon2Gaps[j-1]))
                duplicationSizes, duplicationDetails, operon2Gaps[j-1], operon2GapPositions[j-1] = checkForMatchesWithinOperons(strain2.genomeFragments, event.fragmentDetails2, operon2Gaps[j-1], operon2GapPositions[j-1])
                strain2 = addDuplicationEventsToStrain(strain2, duplicationSizes, duplicationDetails) #Adds duplication details to strain

                #incrementDuplicateSizeCounters(duplicationSizes)
                #incrementDeletionSizeCounters(deletionSizes)
                #print('Gap being processed: %s' % (operon2Gaps[j-1]))
                #print('Number of unique genes found: %s' %(numUniqueFound))
                #print('Number of deletion genes found: %s' %(deletionSizes))
                #print('Number of duplicate genes found: %s' %(duplicationSizes))
                if len(operon2Gaps[j-1]) > 0:
                    deletionDetails = ''
                    deletionSizes = []
                    #Insert gap into operon
                    operon2Gaps[j-1].reverse()
                    for k in range (0, len(operon2Gaps[j-1])):
                        ancestralOperon.insert(operon2GapIndexes[j-1], operon2Gaps[j-1][k])
                        deletionDetails += operon2Gaps[j-1][k] + ' ' + str(operon2GapPositions[j-1][k] + event.fragmentDetails1.startPositionInGenome) + ',' #TODO this might change if operon is in negative orientation
                    deletionDetails += '|'                      #End of deleted segment
                    deletionSizes.append(len(operon2Gaps[j-1])) #Size of segment
                    strain2 = addDeletionEventsToStrain(strain2, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
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
# checkForMatchesWithinOperons
# Parameters:
# Description: Checks if a given gap exists in any of the other operons
######################################################
def checkForMatchesWithinOperons(genomeFragments, fragment, gap, positions):
    allDuplicationSizes = []
    allDuplicationDetails = ''

    for w in range(0, len(genomeFragments)): #Iterate through all fragments
        currFragment = genomeFragments[w]
        if currFragment.startPositionInGenome != fragment.startPositionInGenome and len(gap) > 0 and len(currFragment.sequence) > 1: #They're not the same operon and the gap is longer than 0 and the operon is not a singleton
            duplicationSizes, duplicationDetails, gap, positions = checkForMatch(gap, positions, currFragment.sequence, fragment)

            if len(duplicationSizes) > 0: #If we found a duplicate, add it to the totals
                allDuplicationSizes.append(duplicationSizes)
                allDuplicationDetails += duplicationDetails

    return allDuplicationSizes, allDuplicationDetails, gap, positions

######################################################
# checkForMatch
# Parameters: gap: genes in gap, sequence: sequence to check if the genes exist there
# Description: Checks if a given gap exists in a given sequence
######################################################
def checkForMatch(gap, positions, sequence, fragment):
    geneDuplicateSizes = [] #An array of duplicate sizes
    duplicationDetails = '' #Details about the duplication
    windowSize = len(gap)   #Size of window
    startIndex = 0          #Start position of window
    endIndex = len(gap)     #End position of Window

    while windowSize > 1:
        genes = gap[startIndex:endIndex] #Grabs the genes within the window
        genesMatched = 0

        for x in range(0, len(sequence)): #Iterate over all genes in sequence
            if genes[0] == sequence[x] and genesMatched == 0: #If the first gene matches in the sequence and we don't have an existing match then check if this is potentially a match for the gap
                for y in range(0, len(genes)): #Iterate over all of the genes in the gap
                    if (x+y) < len(sequence) and genes[y] == sequence[y+x]: #Match sure we don't go out of bounds and the genes match
                        genesMatched+=1
                if genesMatched != len(genes): #Reset the counter if we don't find a match for the gap
                    genesMatched = 0

        if genesMatched != 0 and genesMatched == len(genes): #If we found a match for the gap in the sequence
            geneDuplicateSizes.append(len(genes))
            duplicateGenes = gap[startIndex:endIndex]
            duplicatePositions = positions[startIndex:endIndex]

            #Creates a string containing the genes and their position
            for p in range(0, len(duplicateGenes)):
                gene = duplicateGenes[p]
                pos = duplicatePositions[p]
                duplicationDetails += gene + ' ' + str(pos + fragment.startPositionInGenome) + ', ' #TODO might need some special calculation if the operon is reversed!
            duplicationDetails += '|' #This indicates end of duplication fragment

            #Remove the duplicated region
            del gap[startIndex:endIndex]
            del positions[startIndex:endIndex]
            startIndex = endIndex
        else:
            startIndex+=1

        if (startIndex + windowSize) > len(gap):
            #reduce and reset starting position of the window
            windowSize = min(windowSize-1, len(gap))
            startIndex = 0
            endIndex = startIndex + windowSize

    return geneDuplicateSizes, duplicationDetails, gap, positions

######################################################
# checkForMatchesWithinAlignment
# Parameters:
# Description: Takes an array of gaps and an alignment, then checks if the genes in the gap match any of the genes in the alignment, if they do then genes are popped off the gap array
######################################################
def checkForMatchesWithinAlignment(arrayOfGaps, alignedGenes, arrayOfGapPositions, fragment):
    allDuplicationSizes = []
    allDuplicationDetails = ''

    for w in range(0, len(arrayOfGaps)):
        gap = arrayOfGaps[w] #Genes within gap
        positions = arrayOfGapPositions[w] #Positions of genes within gap
        duplicationSizes, duplicationDetails, gap, positions = checkForMatch(gap, positions, alignedGenes, fragment)

        if len(duplicationSizes) > 0: #If we found a duplicate, add it to the totals
            allDuplicationSizes.append(duplicationSizes)
            allDuplicationDetails += duplicationDetails
            arrayOfGaps[w] = gap
            arrayOfGapPositions[w] = positions

    return arrayOfGaps, arrayOfGapPositions, allDuplicationSizes, allDuplicationDetails