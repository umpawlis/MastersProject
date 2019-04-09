from SequenceService import formatAndComputeOperonDifferences
#from SequenceService import incrementDuplicateSizeCounters
#from SequenceService import incrementDeletionSizeCounters
from SequenceService import addDuplicationEventsToStrain
from SequenceService import addDeletionEventsToStrain
from GlobalAlignmentModule import performGlobalAlignment
from SequenceService import reverseSequence
from Event import Event
import globals
import copy

######################################
###Self Global Alignment Functions####
######################################

######################################################
# findOrthologsBySelfGlobalAlignment
# Parameters:
# Description: Performs a self global alignment to identify duplicate operons
######################################################
def findOrthologsBySelfGlobalAlignment(strain, coverageTracker, targetStrain):
    print('Performing self-global alignment on strain: %s' %(strain.getName()))
    lossEvents = []
    duplicationEvents = []
    sequence = strain.getSequence()

    for x in range(0, len(coverageTracker)):
        #Check if marked
        if coverageTracker[x] == False:
            bestScore = 1000 #Make the best score some large numer
            bestEvent = None #Initialize the event
            minDistance = 1000 #Used to track the minimum distance from singleton to operon that has an identical gene

            if len(sequence[x].split(',')) > 1: #This is an operon
                for y in range(0, len(sequence)):
                    #make sure we're not comparing the same operons and that the second operon is NOT  a singleton
                    if x != y and len(sequence[y].split(',')) > 1 and coverageTracker[y] == True:
                        op1 = sequence[x]
                        op2 = sequence[y]

                        #Gene content differences ie which one has more genes, operons 1 and 2, and the number of unique genes between the operons being compared (ie does the operon have any unique genes)
                        geneContentDifferences, operon1, operon2, numUniqueGenes = formatAndComputeOperonDifferences(op1, op2)
                        #Checks if either operons are in the - orientation
                        negativeOrientationOp1 = reverseSequence(op1)
                        negativeOrientationOp2 = reverseSequence(op2)

                        #Tracks whether we reversed the operons
                        operon1Reversed = False
                        operon2Reversed = False

                        #Create event for this comparison
                        event = Event(0)
                        event.setGenome1Operon(operon1)
                        event.setGenome2Operon(operon2)
                        event.setGenome1Name(strain.getName())
                        event.setGenome2Name(strain.getName())
                        event.isOriginallyNegativeOrientationOp1(negativeOrientationOp1) #This tracks the original orientation of op1
                        event.isOriginallyNegativeOrientationOp2(negativeOrientationOp2) #This tracks the original orientation of op2
                        event.setOperon1Index(x)
                        event.setOperon2Index(y)
                        event.setTechnique('Self Global Alignment (Operon)')

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
                        event.setAncestralOperonGeneSequence(copy.deepcopy(operon1)) #Set the ancestral operon sequence
                        score, event = performGlobalAlignment(operon1, operon2, event) #Perform the global alignment
                        event.setScore(score)

                        threshold = max(len(operon1), len(operon2))
                        threshold = threshold//3
                        if geneContentDifferences <= threshold and score < bestScore:
                            bestScore = score
                            bestEvent = event
            else: #This is a singleton gene
                for y in range(0, len(sequence)):
                    if x != y and coverageTracker[y] == True: #Make sure we're not comparing the same singleton genes
                        op1 = sequence[x]
                        op2 = sequence[y]

                        #Gene content differences ie which one has more genes, operons 1 and 2, and the number of unique genes between the operons being compared (ie does the operon have any unique genes)
                        geneContentDifferences, operon1, operon2, numUniqueGenes = formatAndComputeOperonDifferences(op1, op2)
                        #Checks if either operons are in the - orientation
                        negativeOrientationOp1 = reverseSequence(op1)
                        negativeOrientationOp2 = reverseSequence(op2)

                        #Tracks whether we reversed the operons
                        operon1Reversed = False
                        operon2Reversed = False

                        #Create event for this comparison
                        event = Event(0)
                        event.setGenome1Operon(operon1)
                        event.setGenome2Operon(operon2)
                        event.setGenome1Name(strain.getName())
                        event.setGenome2Name(strain.getName())
                        event.isOriginallyNegativeOrientationOp1(negativeOrientationOp1) #This tracks the original orientation of op1
                        event.isOriginallyNegativeOrientationOp2(negativeOrientationOp2) #This tracks the original orientation of op2
                        event.setOperon1Index(x)
                        event.setOperon2Index(y)
                        event.setTechnique('Self Global Alignment (Singleton)')

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
                        event.setAncestralOperonGeneSequence(copy.deepcopy(operon1)) #Set the ancestral operon sequence

                        if operon1[0] in operon2 and abs(x-y) < minDistance: #Checks if the singleton gene is located in the operon and if the distance is smaller
                            minDistance = abs(x-y)
                            event.setScore(0)
                            bestEvent = event

            #Take the event and append it to the duplicate event list
            if bestEvent != None:
                globals.trackingId += 1
                bestEvent.trackingEventId = globals.trackingId
                coverageTracker[x] = True
                duplicationEvents.append(bestEvent)
                
                #Increment the duplicate counter with size of operon since the operon is a duplication
                #incrementDuplicateSizeCounters([len(event.genome1Operon)])
                strain = addDuplicationEventsToStrain([len(event.genome1Operon)], strain)
                
                print('\n&&&&&& Self Global Alignment &&&&&')
                bestEvent.printEvent()
                print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
            else:
                coverageTracker[x] = True
                globals.trackingId +=1
                event = Event(globals.trackingId)
                event.setScore(-1)
                event.setGenome2Operon([])
                event.setGenome1Name(strain.getName())
                event.setGenome2Name('None')
                event.isOriginallyNegativeOrientationOp1(reverseSequence(sequence[x]))
                event.isOriginallyNegativeOrientationOp2(False)
                event.setOperon1Index(x)
                event.setOperon2Index(-1)
                event.setTechnique('Self Global Alignment (No match!)')
                event.isReversedOp1(False)
                event.isReversedOp2(False)

                #Set the ancestral operon sequence
                ancestralGenes = []
                operonGenes = sequence[x]
                operonGenes = operonGenes.replace('-', '')
                operonGenes = operonGenes.replace('[', '')
                operonGenes = operonGenes.replace(']', '')
                operonGenesList = operonGenes.split(',')
                for gene in operonGenesList:
                    ancestralGenes.append(gene.strip())
                event.setAncestralOperonGeneSequence(ancestralGenes)
                event.setGenome1Operon(copy.deepcopy(ancestralGenes))
                lossEvents.append(event)

                #Increment the loss counter with the size of the operon since the operon is a loss
                #incrementDeletionSizeCounters([len(event.genome1Operon)])
                targetStrain = addDeletionEventsToStrain([len(event.genome1Operon)], targetStrain)
                
                print('\n&&&&&& Self Global Alignment &&&&&')
                event.printEvent()
                print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')

    return duplicationEvents, lossEvents, coverageTracker, targetStrain, strain