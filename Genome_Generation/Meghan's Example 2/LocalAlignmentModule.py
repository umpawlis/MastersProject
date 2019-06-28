import globals
import numpy as np
from Event import Event
from GlobalAlignmentModule import constructStatement
from GlobalAlignmentModule import reconstructOperonSequence

#Local alignment parameters
matchWithCodon = 1.0
matchWithoutCodon = 0.5
mismatch = -1.0
gap = -1.0
    
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
    minValue = min(len(coverageTracker1), len(coverageTracker2)) #Used for computing all of the values in the matrix

    #Scan the matrix x times to find optimal local alignments everytime an entire matrix is computed
    for x in range(0, minValue):
        highestScore = -1
        distance = 50 #An arbitrary large number

        for i in range(0, len(coverageTracker1)):
            filteredList = iter(filter(lambda x:x.fragmentIndex == i, strain1.genomeFragments)) #Get the fragment we need based on the index
            fragment1 = next(filteredList, None)

            if coverageTracker1[i] == False and len(fragment1.sequence) > 1: #Make sure the operon is not marked and not a singleton

                for j in range(0, len(coverageTracker2)):
                    filteredList = iter(filter(lambda x:x.fragmentIndex == j, strain2.genomeFragments)) #Get the fragment we need based on the index
                    fragment2 = next(filteredList, None)

                    if coverageTracker2[j] == False and len(fragment2.sequence) > 1: #Make sure the operon is not marked and not a singleton
                        #Initialize the event
                        event = Event(0)
                        event.setFragmentDetails1(fragment1)
                        event.setFragmentDetails2(fragment2)
                        event.setGenome1Name(strain1.name)
                        event.setGenome2Name(strain2.name)
                        event.setTechnique('Local Alignment')
                        
                        #Get the newest copy of the fragment into event as it may have been update
                        filteredList = iter(filter(lambda x : x.fragmentIndex == event.fragmentDetails1.fragmentIndex, strain1.genomeFragments))
                        fragment = next(filteredList, None)
                        event.fragmentDetails1 = fragment
                        filteredList = iter(filter(lambda x : x.fragmentIndex == event.fragmentDetails2.fragmentIndex, strain2.genomeFragments))
                        fragment = next(filteredList, None)
                        event.fragmentDetails2 = fragment
                        
                        score, startPosition, endPosition, event = localAlignment(event.fragmentDetails1, event.fragmentDetails2, event)

                        if (score > highestScore) or (score == highestScore and (abs(i - j)) < distance):
                            highestScore = score
                            rowIndex = i
                            colIndex = j
                            distance = abs(i - j)
                            chosenOpEvent = event
        #After scanning the whole matrix if we found a best score, then store it
        if highestScore > -1:
            print('\n******** Local Alignment ***********')
            print('**************************************\n')
            event = chosenOpEvent
            globals.trackingId += 1
            localAlignmentCounter+=1
            coverageTracker1[rowIndex] = True
            coverageTracker2[colIndex] = True
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

            events.append(event) #Add the event to the tracking events list
            print(event.toString())
            print('\n**************************************')
            print('**************************************\n\n')

    return events, coverageTracker1, coverageTracker2, localAlignmentCounter, strain1, strain2

######################################################
# localAlignment
# Parameters:
# Description: Performs local alignment on two operons
######################################################
def localAlignment(fragment1, fragment2, event):
    operon1 = fragment1.sequence
    operon2 = fragment2.sequence

    #Initialize the matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))
    maxScore = 0
    maxPosition = (0, 0)
    returningScore = -1
    shortestLength = min(len(operon1), len(operon2)) #Used to determine if this comparison is of interest to us

    #Perform the Local Alignment
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            if operon1[a-1].split('_')[0].strip() == operon2[b-1].split('_')[0].strip(): #Genes are identical (without codons)
                if operon1[a-1].strip() == operon2[b-1].strip(): #Genes are identical (with codon)
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + matchWithCodon
                else: #Check all possible options, it was observed that codon mismatches where not always the optimal choice
                    scoreMatrix[a][b] = max(scoreMatrix[a-1][b-1] + matchWithoutCodon, scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)
            else: #Genes don't match
                scoreMatrix[a][b] = max(scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)

            if scoreMatrix[a][b] > maxScore: #Track the position of the max score (that will be our starting point during the trace back)
                maxScore = scoreMatrix[a][b]
                maxPosition = (a, b)
                
    event, endPosition = traceback(operon1, operon2, scoreMatrix, maxPosition, event)
    event.setScore(maxScore)
    
    val1 = len(event.operon1Alignment)
    val2 = (0.75 * shortestLength)
    if val1 >= val2:
        returningScore = maxScore

    return returningScore, maxPosition, endPosition, event

######################################################
# traceback
# Parameters: operon1, operon2, scoreMatrix, startPosition
# Description: treaverses score matrix to determine optimal alignment
######################################################
def traceback(operon1, operon2, scoreMatrix, startPosition, event):
    match = 1
    codonMismatch = 0
    mismatch = 0
    numGaps = 0

    END, DIAG, UP, LEFT = range(4)
    x, y = startPosition
    
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

    #Tracks where the extra genes are from
    gap1Indexes = []
    gap2Indexes = []

    #Tracks the aligned genes
    aligned_seq1 = []
    aligned_seq2 = []

    #Need to capture the genes outside the aligned regions. Remember this traceback only traverses the regions with a score above 0! ie we'll be missing genes
    currIndex = len(operon1)
    while currIndex > 1 and currIndex > x:
        operon2Gap.insert(0, operon1[currIndex - 1])
        operon2GapIndex.insert(0, currIndex - 1)
        currIndex -= 1
    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2GapIndexes.insert(0, operon2GapIndex)
        gap2Indexes.insert(0, len(aligned_seq2))
        operon2Gap = []
        operon2GapIndex = []
    
    currIndex = len(operon2)
    while currIndex > 1 and currIndex > y:
        operon1Gap.insert(0, operon2[currIndex - 1])
        operon1GapIndex.insert(0, currIndex - 1)
        currIndex -= 1
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1GapIndexes.insert(0, operon1GapIndex)
        gap1Indexes.insert(0, len(aligned_seq1))
        operon1Gap = []
        operon1GapIndex = []
    
    move = nextMove(scoreMatrix, x, y, operon1[x-1], operon2[y-1])          
    while move != END:
        if move == DIAG:
            if operon1[x-1] == operon2[y-1]: #Both the gene and codon match
                match += 1
            elif operon1[x-1].split('_')[0].strip() == operon2[y-1].split('_')[0].strip(): #Codon mismatch
                codonMismatch += 1
                codonMismatchIndexesStrain1.append(x-1)
                codonMismatchGenesStrain1.append(operon1[x-1])
                codonMismatchIndexesStrain2.append(y-1)
                codonMismatchGenesStrain2.append(operon2[y-1])
            else: #Substitution
                mismatch += 1
                substitutionIndexesStrain1.append(x-1)
                substitutionGenesStrain1.append(operon1[x-1])
                substitutionIndexesStrain2.append(y-1)
                substitutionGenesStrain2.append(operon2[y-1])

            aligned_seq1.insert(0, operon1[x - 1]) #insert gene into alignment
            aligned_seq2.insert(0, operon2[y - 1]) #insert gene into alignment
            x -= 1 #move x to the next position
            y -= 1 #move y to the next position
            operon1ConsecutiveGap = False #Reset the consecutive gap tracker for operon 1
            operon2ConsecutiveGap = False #Reset the consecutive gap tracker for operon 2

        elif move == UP:
            index = x - 1
            # aligned_seq1.insert(0, operon1[x - 1])
            # aligned_seq2.insert(0, ' - ')
            numGaps += 1
            x -= 1

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
                gap2Indexes.insert(0, len(aligned_seq2))
                operon2ConsecutiveGap = True

        else:
            index = y - 1
            # aligned_seq1.insert(0, ' - ')
            # aligned_seq2.insert(0, operon2[y - 1])
            numGaps += 1
            y -= 1

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
                gap1Indexes.insert(0, len(aligned_seq1))
                operon1ConsecutiveGap = True

        move = nextMove(scoreMatrix, x, y, operon1[x-1], operon2[y-1]) #Makes move to the next cell in the matrix

    aligned_seq1.insert(0, operon1[x - 1]) #Add the last aligned gene
    aligned_seq2.insert(0, operon2[y - 1])
    endPosition = (x,y)

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
        
    #Need to capture the genes outside the aligned regions. Remember this traceback only traverses the regions with a score above 0! ie we'll be missing genes
    currIndex = x - 1
    while currIndex > 0:
        operon2Gap.insert(0, operon1[currIndex - 1])
        operon2GapIndex.insert(0, currIndex - 1)
        currIndex -= 1
    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2GapIndexes.insert(0, operon2GapIndex)
        gap2Indexes.insert(0, len(aligned_seq2))
        operon2Gap = []
        operon2GapIndex = []
    
    currIndex = y - 1
    while currIndex > 0:
        operon1Gap.insert(0, operon2[currIndex - 1])
        operon1GapIndex.insert(0, currIndex - 1)
        currIndex -= 1
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1GapIndexes.insert(0, operon1GapIndex)
        gap1Indexes.insert(0, len(aligned_seq1))
        operon1Gap = []
        operon1GapIndex = []
                
    #The indexes values need to be flipped b/c right now they're oriented from right to left
    if len(gap1Indexes) > 0:
        for x in range(0, len(gap1Indexes)):
            gap1Indexes[x] = len(aligned_seq1) - gap1Indexes[x]
    if len(gap2Indexes) > 0:
        for x in range(0, len(gap2Indexes)):
            gap2Indexes[x] = len(aligned_seq2) - gap2Indexes[x]

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
    event.setNumMismatches(numGaps)

    #Codon Mismatch details
    event.setNumCodonMismatches(codonMismatch)
    event.setCodonMismatchGenesStrain1(codonMismatchGenesStrain1)
    event.setCodonMismatchGenesStrain2(codonMismatchGenesStrain2)
    event.setCodonMismatchIndexesStrain1(codonMismatchIndexesStrain1)
    event.setCodonMismatchIndexesStrain2(codonMismatchIndexesStrain2)

    #Substitution details
    event.setNumSubstitutions(mismatch)
    event.setSubstitutionIndexesStrain1(substitutionIndexesStrain1)
    event.setSubstitutionIndexesStrain2(substitutionIndexesStrain2)
    event.setSubstitutionGenesStrain1(substitutionGenesStrain1)
    event.setSubstitutionGenesStrain2(substitutionGenesStrain2)

    #Alignment details
    event.setOperon1Alignment(aligned_seq1)
    event.setOperon2Alignment(aligned_seq2)

    #Gap details
    event.setOperon1Gaps(operon1Gaps)
    event.setOperon2Gaps(operon2Gaps)
    event.setOperon1GapPositions(operon1GapIndexes)
    event.setOperon2GapPositions(operon2GapIndexes)
    event.setOperon1GapIndexes(gap1Indexes)
    event.setOperon2GapIndexes(gap2Indexes)

    return event, endPosition

######################################################
# nextMove
# Parameters: scoreMatrix, x, y
# Description: determines which direction to traverse score matrix
######################################################
def nextMove(scoreMatrix, x, y, gene1, gene2):
    diag = scoreMatrix[x - 1][y - 1]
    up   = scoreMatrix[x - 1][y]
    left = scoreMatrix[x][y - 1]
    
    gene1Array = gene1.strip().split('_')
    gene2Array = gene2.strip().split('_')
    
    if (diag + matchWithCodon) == scoreMatrix[x][y] and gene1.strip() == gene2.strip():     #Perfect match
        return 1 if diag !=0 else 0
    elif (diag + matchWithoutCodon) == scoreMatrix[x][y] and len(gene1Array) == 2 and len(gene2Array) == 2 and gene1Array[0] == gene2Array[0] and gene1Array[1] != gene2Array[1]: #Codon mismatch
        return 1 if diag !=0 else 0
    elif (diag + mismatch) == scoreMatrix[x][y] and gene1Array[0] != gene2Array[0]:         #Substitution
        return 1 if diag !=0 else 0
    elif (left + gap) == scoreMatrix[x][y]:                                                 #Gap
        return 3 if left != 0 else 0
    elif (up + gap) == scoreMatrix[x][y]:                                                   #Gap
        return 2 if up != 0 else 0
    else:
        return 0