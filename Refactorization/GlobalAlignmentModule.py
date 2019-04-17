from Event import Event
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
    #events, coverageTracker1, coverageTracker2, globalAlignmentCounter, strain1, strain2 = scanGlobalAlignmentMatrixForOrthologs(alignmentMatrix, eventMatrix, coverageTracker1, coverageTracker2, strain1, strain2)

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

            #Case 1: An origin or terminus
            if (op1.description == 'Origin' and op2.description == 'Origin') or (op1.description == 'Terminus' and op2.description == 'Terminus'):

                event = Event(0)
                event.setScore(0.0)
                event.setFragmentDetails1(op1)
                event.setFragmentDetails2(op2)
                event.setGenome1Name(strain1.name)
                event.setGenome2Name(strain2.name)
                event.setTechnique('Global Alignment')
                event.setOperon1Alignment(copy.deepcopy(op1.sequence))
                event.setOperon2Alignment(copy.deepcopy(op2.sequence))
                
                eventMatrix[x][y] = event
                
                globalAlignmentMatrix[x][y] = str(0) + '*'

            #Case 2: Two singleton genes are being compared
            elif len(op1.sequence) == 1 and len(op2.sequence) == 1:

                if op1.sequence[0] == op2.sequence[0]:
                    event = Event(0)
                    event.setScore(0.0)
                    event.setFragmentDetails1(op1)
                    event.setFragmentDetails2(op2)
                    
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
                    
                    eventMatrix[x][y] = event
                    globalAlignmentMatrix[x][y] = str(0) + '*' #Singletons are a perfect match
                else:
                    globalAlignmentMatrix[x][y] = -999 #Singletons don't match

            #Case 3: Both are operons
            elif len(op1.sequence) > 1 and len(op2.sequence) > 1:
                event = Event(0)
                event.setFragmentDetails1(op1)
                event.setFragmentDetails2(op2)
                
                event.setGenome1Operon(op1.sequence)
                event.setGenome2Operon(op2.sequence)
                event.setGenome1Name(strain1Name)
                event.setGenome2Name(strain2Name)
                event.isOriginallyNegativeOrientationOp1(op1.isNegativeOrientation)
                event.isOriginallyNegativeOrientationOp2(op2.isNegativeOrientation)
                event.setOperon1Index(x)
                event.setOperon2Index(y)
                event.setTechnique('Global Alignment')
                

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