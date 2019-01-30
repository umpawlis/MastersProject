from SequenceService import formatAndComputeOperonDifferences
from SequenceService import reverseSequence
from Event import Event

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
    
    #trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter = scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, operonEventMatrix, coverageTracker1, coverageTracker2, strain1, strain2)

    return events, coverageTracker1, coverageTracker2, globalAlignmentCounter

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

            #Gene content differences ie which one has more genes, operons 1 and 2, and the number of unique genes between the operons being compared (ie does the operon have any unique genes)
            geneContentDifferences, operon1, operon2, numUniqueGenes = formatAndComputeOperonDifferences(op1, op2)

            #Checks if either operons are in the - orientation
            negativeOrientationOp1 = reverseSequence(op1)
            negativeOrientationOp2 = reverseSequence(op2)
            
            #Tracks whether we reversed the operons
            operon1Reversed = False
            operon2Reversed = False

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
                
            #Case 1: We have two singleton operons
            if len(operon1) == 1 and len(operon2) == 1:
                #Perfect match
                if operon1 == operon2:
                    globalAlignmentMatrix[x][y] =  str(0) + '*'
                    
                    event = Event(0)
                    event.setScore(0)
                    
                    
                    
                    #We need this information to track whether the operon was in the negative orientation and if it was flipped
                    operonEvents = OperonEvents(1, 0, 0, 0, operon1, operon2, None)
                    operonEvents.setOperon1Reversed(operon1Reversed)
                    operonEvents.setOperon2Reversed(operon2Reversed)
                    operonEvents.setOperon1NegativeOrientation(operon1NegativeOrientation)
                    operonEvents.setOperon2NegativeOrientation(operon2NegativeOrientation)
                    
                    operonEventMatrix[x][y] = operonEvents

                #Mismatch
                else:
                    globalAlignmentMatrix[x][y] = -999

            #Case 2: Only one of them is a singleton operon
            elif (len(operon1) == 1 and len(operon2) > 1) or (len(operon2) == 1 and len(operon1) > 1):
                globalAlignmentMatrix[x][y] = -999

            #Case 3: None of them are singleton operons, perform a global alignment
            elif len(op1) > 1 and len(op2) > 1:
                score, operonEvents = performGlobalAlignment(operon1, operon2)

                #We need this information to track whether the operon was in the negative orientation and if it was flipped
                operonEvents.setOperon1Reversed(operon1Reversed)
                operonEvents.setOperon2Reversed(operon2Reversed)
                operonEvents.setOperon1NegativeOrientation(operon1NegativeOrientation)
                operonEvents.setOperon2NegativeOrientation(operon2NegativeOrientation)

                globalAlignmentMatrix[x][y] = score
                operonEventMatrix[x][y] = operonEvents

                threshold = max(len(operon1), len(operon2))
                threshold = threshold//3

                if geneContentDifferences <= threshold:
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