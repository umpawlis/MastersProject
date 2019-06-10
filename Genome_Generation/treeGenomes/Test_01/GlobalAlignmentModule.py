from SequenceService import addDuplicationEventsToStrain
from SequenceService import addDeletionEventsToStrain
from SequenceService import computeOperonDifferences
from DeletionDetails import DeletionDetails
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
            event.setScore(1.0)
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
                globalAlignmentMatrix[x][y] = str(1) + '*'

            #Case 2: Two singleton genes are being compared and not the origin or terminus
            elif len(op1.sequence) == 1 and len(op2.sequence) == 1 and op1.description != 'Origin' and op2.description != 'Origin' and op1.description != 'Terminus' and op2.description != 'Terminus':
                if op1.sequence[0] == op2.sequence[0]:
                    event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                    event.setOperon2Alignment(copy.deepcopy(op2.sequence))

                    eventMatrix[x][y] = event
                    globalAlignmentMatrix[x][y] = str(1) + '*'
                else:
                    globalAlignmentMatrix[x][y] = -999 #Singletons don't match

            #Case 3: Both are operons
            elif len(op1.sequence) > 1 and len(op2.sequence) > 1:
                redoAlignment1 = False #Tracks whether we should redo the alignment to potentially get a better score
                redoAlignment2 = False #Tracks whether we should redo the alignment to potentially get a better score
                firstAttemptEvent = copy.deepcopy(event)
                firstAttemptScore, firstAttemptEvent = performGlobalAlignment(op1.sequence, op2.sequence, firstAttemptEvent)
                
                #Remove genes that were deleted twice in a row (check both fragments)
                if len(event.fragmentDetails1.deletionDetailsList) > 0:
                    eventCopy1 = copy.deepcopy(firstAttemptEvent)
                    op1Copy = copy.deepcopy(op1)
                    redoAlignment1, op1Copy, deletionDetailsList1 = removeGenesDeletedMultipleGenerations(eventCopy1, op1Copy, eventCopy1.operon1Gaps, eventCopy1.operon1GapPositions, eventCopy1.fragmentDetails1.fragmentIndex)
                if len(event.fragmentDetails2.deletionDetailsList) > 0:
                    eventCopy2 = copy.deepcopy(firstAttemptEvent)
                    op2Copy = copy.deepcopy(op2)
                    redoAlignment2, op2Copy, deletionDetailsList2 = removeGenesDeletedMultipleGenerations(eventCopy2, op2Copy, eventCopy2.operon2Gaps, eventCopy2.operon2GapPositions, eventCopy2.fragmentDetails2.fragmentIndex)
                
                #Redo the Alignment if any genes were removed from the operon
                if redoAlignment1 or redoAlignment2:
                    secondAttemptEvent = copy.deepcopy(event)
                    if redoAlignment1:
                        secondAttemptOp1 = op1Copy
                        secondAttemptEvent.setFragmentDetails1(op1Copy)
                        secondAttemptEvent.fragmentDetails1.deletionDetailsList = deletionDetailsList1
                    else:
                        secondAttemptOp1 = op1
                    if redoAlignment2:
                        secondAttemptOp2 = op2Copy
                        secondAttemptEvent.setFragmentDetails2(op2Copy)
                        secondAttemptEvent.fragmentDetails2.deletionDetailsList = deletionDetailsList2
                    else:
                        secondAttemptOp2 = op2
                    secondAttemptScore, secondAttemptEvent = performGlobalAlignment(secondAttemptOp1.sequence, secondAttemptOp2.sequence, secondAttemptEvent)
                    
                    if secondAttemptScore > firstAttemptScore:
                        bestScore = secondAttemptScore
                        secondAttemptEvent.genesDeletedFromOperon = True
                        bestEvent = secondAttemptEvent
                    else:
                        bestScore = firstAttemptScore
                        bestEvent = firstAttemptEvent
                else:
                    bestScore = firstAttemptScore
                    bestEvent = firstAttemptEvent
                #At this point we have the best scores
                bestEvent.setScore(bestScore)
                eventMatrix[x][y] = bestEvent
                globalAlignmentMatrix[x][y] = bestScore

                #threshold = max(len(op1.sequence), len(op2.sequence))
                #threshold = threshold//3
                #numOperonDifferences = computeOperonDifferences(op1.sequence, op2.sequence)
                #if numOperonDifferences <= threshold:
                if bestScore >= 0:
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
# checkDeletionEvents
# Parameters:
# Description: Checks if given gene and position exists in the deletion details list
######################################################
def removeGenesDeletedMultipleGenerations(event, op, operonGaps, operonGapPositions, fragmentId):
    if len(operonGaps) > 0:
        operonGaps.reverse()             #Reverse the arrays so we start removing from the highest positions
        operonGapPositions.reverse()     #Reverse the arrays so we start removing from the highest positions
        sequenceChanged = False          #Tracks whether we actually removed any genes
        for t in range(0, len(operonGaps)):
            genes = operonGaps[t]
            positions = operonGapPositions[t]
            for g in range(0, len(genes)):
                gene = genes[g]
                position = positions[g]
                #Check if this combination exists in the deletion tracker
                removeGenes = checkDeletionEvents(gene, position, event.fragmentDetails1.deletionDetailsList, fragmentId)
                if removeGenes == True:
                    sequenceChanged = True
                    redoAlignment = True
                    del op.sequence[position]
        #Reconstruct the operon string sequence now
        if sequenceChanged:
            sequenceCopy = copy.deepcopy(op.sequence)
            stringOperon = ''
            if op.isNegativeOrientation == True:
                sequenceCopy.reverse()
                stringOperon += '-'
            stringOperon += '['
            for q in range(0, len(sequenceCopy)):
                geneCopy = sequenceCopy[q]
                if q != (len(sequenceCopy) - 1):
                    stringOperon += geneCopy + ', '
                else:
                    stringOperon += geneCopy
            stringOperon += ']'
            op.originalSequence = stringOperon
    return redoAlignment, op, event.fragmentDetails1.deletionDetailsList

######################################################
# checkDeletionEvents
# Parameters:
# Description: Checks if given gene and position exists in the deletion details list
######################################################
def checkDeletionEvents(gene, position, deletionList, fragmentId):
    for x in range(0, len(deletionList)):
        deletion = deletionList[x]
        #Make sure these are the genes from the same operon!!
        if deletion.ancestralGene == gene and int(deletion.ancestralPosition) == int(position) and fragmentId == deletion.ancestralFragmentId:
            deletion.geneRemoved = True
            return True
    return False

######################################################
# performGlobalAlignment
# Parameters:
# Description: Performs a global alignment
######################################################
def performGlobalAlignment(operon1, operon2, event):

    #initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))

    for a in range(0, len(operon1)+1):
        scoreMatrix[a][0] = -a

    for a in range(0, len(operon2)+1):
        scoreMatrix[0][a] = -a

    #perform the Global Alignment
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            #check if genes are identical
            if operon1[a-1].split('_')[0].strip() == operon2[b-1].split('_')[0].strip():
                #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                if operon1[a-1].strip() == operon2[b-1].strip():
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + globals.match
                else:
                    #Solves a special case with a bunch of duplicates with different codons
                    scoreMatrix[a][b] = max(scoreMatrix[a-1][b-1] + globals.codonCost, scoreMatrix[a-1][b] + globals.deletionCost, scoreMatrix[a][b-1] + globals.deletionCost, scoreMatrix[a-1][b-1] + globals.substitutionCost)
            else:
                scoreMatrix[a][b] = max(scoreMatrix[a-1][b] + globals.deletionCost, scoreMatrix[a][b-1] + globals.deletionCost, scoreMatrix[a-1][b-1] + globals.substitutionCost)

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
        if i > 0 and j > 0 and matrix[i][j] == (matrix[i-1][j-1] + globals.match) and operon1[i-1] == operon2[j-1]:
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
    currentScoreSelected = maxValue
    #Keep iterating util we find all the optimal scores (Finding orthologs using global alignment)
    while currentScoreSelected >= 0:
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
                        
                        #Check if we have to go back and remove genes
                        if event.genesDeletedFromOperon == True:
                            #Check which of the strains we need to do a backtrack on
                            doBackTrack1, numGenesDeleted1 = operonHadGenesRemoved(event.fragmentDetails1.deletionDetailsList)
                            doBackTrack2, numGenesDeleted2 = operonHadGenesRemoved(event.fragmentDetails2.deletionDetailsList)
                            if doBackTrack1:
                                #Get strain 1 from global strains array
                                filteredList = iter(filter(lambda x : x.name == event.genome1Name, globals.strains))
                                ancestor1 = next(filteredList, None)
                                if ancestor1 != None:
                                    #Get the operon that was modified
                                    filteredList = iter(filter(lambda x : x.fragmentIndex >= event.fragmentDetails1.fragmentIndex, ancestor1.genomeFragments))
                                    fragment = next(filteredList, None)
                                    while fragment:
                                        if fragment.fragmentIndex == event.fragmentDetails1.fragmentIndex:
                                            #Replace the operon array and the string sequence
                                            fragment.originalSequence = copy.deepcopy(event.fragmentDetails1.originalSequence)
                                            fragment.sequence = copy.deepcopy(event.fragmentDetails1.sequence)
                                        else:
                                            #Update the start position of the other operons by subtracting by the number of genes deleted
                                            fragment.startPositionInGenome = int(fragment.startPositionInGenome) - numGenesDeleted1
                                        fragment = next(filteredList, None)
                                    #Now remove one gene at a time from the strains
                                    
                                else:
                                    print('Error! Something went wrong because we failed to find ancestor %s!' % (event.genome1Name))

                        
                        event.trackingEventId = globals.trackingId
                        event, strain1, strain2 = reconstructOperonSequence(event, strain1, strain2)

                        #Codon mismatches
                        if event.numCodonMismatches > 0:
                            codonMismatchDescription1 = constructStatement(event.codonMismatchIndexesStrain1, event.codonMismatchGenesStrain1, event.fragmentDetails1)
                            codonMismatchDescription2 = constructStatement(event.codonMismatchIndexesStrain2, event.codonMismatchGenesStrain2, event.fragmentDetails2)

                            strain1.addCodonMismatchDetails(codonMismatchDescription1)
                            strain2.addCodonMismatchDetails(codonMismatchDescription2)

                        #Substitutions
                        if event.numSubstitutions > 0:
                            substitutionDescription1 = constructStatement(event.substitutionIndexesStrain1, event.substitutionGenesStrain1, event.fragmentDetails1)
                            substitutionDescription2 = constructStatement(event.substitutionIndexesStrain2, event.substitutionGenesStrain2, event.fragmentDetails2)

                            strain1.addSubstitutionDetails(substitutionDescription1)
                            strain2.addSubstitutionDetails(substitutionDescription2)

                        print(event.toString())

                        #Add the event to the tracking events list
                        events.append(event)
                        print('###################################\n')

        currentScoreSelected += (-0.5)
    return events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2

######################################################
# operonHadGenesRemoved
# Parameters:
# Description: Checks an operons deleted gene list if any of the genes were switch from deletions to duplications
######################################################
def operonHadGenesRemoved(deletions):
    count = 0
    result = False
    if len(deletions) > 0:
        for deletion in deletions:
            if deletion.geneRemoved == True:
                result = True
                count += 1
    return result, count

######################################################
# constructStatement
# Parameters:
# Description: constructs the details of where the codon mismatch or substitution occured in the geneome
######################################################
def constructStatement(indexes, genes, fragmentDetails):
    temp = ''
    
    if len(indexes) != len(genes):
        print('Error! These two arrays should be the same length for Codon Mismatch Substitution parallel arrays')
    else:
        for x in range(0, len(indexes)):
            position = indexes[x]
            gene = genes[x]
            
            if fragmentDetails.isNegativeOrientation == False: #Computes position based on whether the operon was originally in the negative orientation
                genePos = position + fragmentDetails.startPositionInGenome
            else:
                genePos = fragmentDetails.startPositionInGenome + len(fragmentDetails.sequence) - position - 1
            
            temp+=gene + ' ' + str(genePos) + ';'

    return temp

######################################################
# reconstructOperonSequence
# Parameters:
# Description: Reconstructs the ancestral operon by determining whether the gaps are losses or duplications
######################################################
def reconstructOperonSequence(event, strain1, strain2):
    ancestralOperon = copy.deepcopy(event.operon1Alignment) #The alignments will be identical except when codon mismatches or substitutions occur
    
    operon1Gaps = event.operon1Gaps
    operon2Gaps = event.operon2Gaps
    
    if len(operon1Gaps) == 0 and len(operon2Gaps) == 0:
        print('No differences detected between these two operons')
        event.setAncestralOperonGeneSequence(ancestralOperon)
    else:
        print('Differences detected between these two operons!')
        
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
        strain2 = addDuplicationEventsToStrain(strain2, duplicateSizesWithinAlignment2, duplicationDetails2)

        #Step 2: Check if these extra genes are the result of a duplication event within another operon, remove them if they are, else insert them
        i = len(operon1Gaps)
        j = len(operon2Gaps)
        
        #Tracks the genes marked as deleted and their positions
        deletedGenes = []
        deletedGenesPositions = []
        fragmentIds = []
        strains = []
        originalDeletedGenes = []
        originalDeletedGenesPositions = []

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
                        
                        if event.fragmentDetails1.isNegativeOrientation == False: #This compute the correct gene position based on whether operon was in the negative orientation or not originally
                            genePos = operon1GapPositions[i-1][k] + event.fragmentDetails1.startPositionInGenome
                        else:
                            genePos = event.fragmentDetails1.startPositionInGenome + len(event.fragmentDetails1.sequence) - operon1GapPositions[i-1][k] - 1
                        
                        #Deleted genes in the ancestral operon
                        deletedGenes.append(operon1Gaps[i-1][k])
                        deletedGenesPositions.append(len(ancestralOperon) - int(operon1GapIndexes[i-1]))
                        fragmentIds.append(event.fragmentDetails1.fragmentIndex)
                        strains.append(strain1)
                        originalDeletedGenes.append(operon1Gaps[i-1][k])
                        originalDeletedGenesPositions.append(genePos)
                        
                        deletionDetails += operon1Gaps[i-1][k] + ' ' + str(genePos) + ', '
                    deletionDetails = deletionDetails[0:(len(deletionDetails) - 2)]
                    deletionDetails += ';'                      #End of deleted segment
                    deletionSizes.append(len(operon1Gaps[i-1])) #Size of segment
                    strain2 = addDeletionEventsToStrain(strain2, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
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
                        
                        if event.fragmentDetails2.isNegativeOrientation == False: #This compute the correct gene position based on whether operon was in the negative orientation or not originally
                            genePos = operon2GapPositions[j-1][k] + event.fragmentDetails2.startPositionInGenome
                        else:
                            genePos = event.fragmentDetails2.startPositionInGenome + len(event.fragmentDetails2.sequence) - operon2GapPositions[j-1][k] - 1
                        
                        #Deleted genes in the ancestral operon
                        deletedGenes.append(operon2Gaps[j-1][k])
                        deletedGenesPositions.append(len(ancestralOperon) - int(operon2GapIndexes[j-1]))
                        fragmentIds.append(event.fragmentDetails2.fragmentIndex)
                        strains.append(strain2)
                        originalDeletedGenes.append(operon2Gaps[j-1][k])
                        originalDeletedGenesPositions.append(genePos)
                        
                        deletionDetails += operon2Gaps[j-1][k] + ' ' + str(genePos) + ', '
                    deletionDetails = deletionDetails[0:(len(deletionDetails) - 2)]
                    deletionDetails += ';'                      #End of deleted segment
                    deletionSizes.append(len(operon2Gaps[j-1])) #Size of segment
                    strain1 = addDeletionEventsToStrain(strain1, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
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
                        
                        if event.fragmentDetails1.isNegativeOrientation == False: #This compute the correct gene position based on whether operon was in the negative orientation or not originally
                            genePos = operon1GapPositions[i-1][k] + event.fragmentDetails1.startPositionInGenome
                        else:
                            genePos = event.fragmentDetails1.startPositionInGenome + len(event.fragmentDetails1.sequence) - operon1GapPositions[i-1][k] - 1
                        
                        #Deleted genes in the ancestral operon
                        deletedGenes.append(operon1Gaps[i-1][k])
                        deletedGenesPositions.append(len(ancestralOperon) - int(operon1GapIndexes[i-1]))
                        fragmentIds.append(event.fragmentDetails1.fragmentIndex)
                        strains.append(strain1)
                        originalDeletedGenes.append(operon1Gaps[i-1][k])
                        originalDeletedGenesPositions.append(genePos)
                        
                        deletionDetails += operon1Gaps[i-1][k] + ' ' + str(genePos) + ', '
                    deletionDetails = deletionDetails[0:(len(deletionDetails) - 2)]
                    deletionDetails += ';'                      #End of deleted segment
                    deletionSizes.append(len(operon1Gaps[i-1])) #Size of segment
                    strain2 = addDeletionEventsToStrain(strain2, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
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
                        
                        if event.fragmentDetails2.isNegativeOrientation == False: #This compute the correct gene position based on whether operon was in the negative orientation or not originally
                            genePos = operon2GapPositions[j-1][k] + event.fragmentDetails2.startPositionInGenome
                        else:
                            genePos = event.fragmentDetails2.startPositionInGenome + len(event.fragmentDetails2.sequence) - operon2GapPositions[j-1][k] - 1
                        
                        #Deleted genes in the ancestral operon
                        deletedGenes.append(operon2Gaps[j-1][k])
                        deletedGenesPositions.append(len(ancestralOperon) - int(operon2GapIndexes[j-1]))
                        fragmentIds.append(event.fragmentDetails2.fragmentIndex)
                        strains.append(strain2)
                        originalDeletedGenes.append(operon2Gaps[j-1][k])
                        originalDeletedGenesPositions.append(genePos)
                        
                        deletionDetails += operon2Gaps[j-1][k] + ' ' + str(genePos) + ', '
                    deletionDetails = deletionDetails[0:(len(deletionDetails) - 2)]
                    deletionDetails += ';'                      #End of deleted segment
                    deletionSizes.append(len(operon2Gaps[j-1])) #Size of segment
                    strain1 = addDeletionEventsToStrain(strain1, deletionSizes, deletionDetails) #Remember, if the genes are detected a deletions, it means it was lost in the other strain!!
                j = j - 1
        
        #Adds deletion details to a global tracker
        if len(deletedGenes) > 0:
            for x in range(0, len(deletedGenes)):
                gene = deletedGenes[x] #Gene that was deleted
                position = len(ancestralOperon) - deletedGenesPositions[x]  #Postion of the gene that was deleted with respect to ancestral operon
                fragmentId = fragmentIds[x]
                strain = strains[x]
                originalGene = originalDeletedGenes[x]
                originalPosition = originalDeletedGenesPositions[x]
                #Construct the deletion tracker object
                details = DeletionDetails(gene, position, fragmentId, strain, originalGene, originalPosition, ancestralOperon)
                event.deletionDetailsList.append(details)
                
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
                allDuplicationSizes.extend(duplicationSizes)
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
                
                if fragment.isNegativeOrientation == False: #Computes position based on whether the operon was originally in the negative orientation
                    genePos = pos + fragment.startPositionInGenome
                else:
                    genePos = fragment.startPositionInGenome + len(fragment.sequence) - pos - 1
                
                duplicationDetails += gene + ' ' + str(genePos) + ', '
            duplicationDetails = duplicationDetails[0:(len(duplicationDetails) - 2)]
            duplicationDetails += ';' #This indicates end of duplication fragment

            #Remove the duplicated region
            del gap[startIndex:endIndex]
            del positions[startIndex:endIndex]
            startIndex = endIndex
            endIndex += windowSize
            #For making sure everything is working
            print('Gap sequence found! Removing gap sequence!')
            print('This is the gap that was removed, %s' % (duplicationDetails))
            print('This is the sequence the gap was found in, %s' % (sequence))
        else:
            startIndex+=1
            endIndex+=1

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
            allDuplicationSizes.extend(duplicationSizes)
            allDuplicationDetails += duplicationDetails
            arrayOfGaps[w] = gap
            arrayOfGapPositions[w] = positions

    return arrayOfGaps, arrayOfGapPositions, allDuplicationSizes, allDuplicationDetails