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
                        #TODO deal with ancestral reconstruction
                        #event, strain1, strain2 = reconstructOperonSequence(event, strain1, strain2)

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