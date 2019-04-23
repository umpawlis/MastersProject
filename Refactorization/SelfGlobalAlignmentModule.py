from SequenceService import addDuplicationEventsToStrain
from GlobalAlignmentModule import performGlobalAlignment
from SequenceService import addDeletionEventsToStrain
from SequenceService import computeOperonDifferences
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
def findOrthologsBySelfGlobalAlignment(strain, coverageTracker):
    print('Performing self-global alignment on strain: %s' %(strain.name))
    lossEvents = []
    duplicationEvents = []
    fragments = strain.genomeFragments
    
    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False: #Operon has not been marked
            bestScore = 1000    #Make the best score some large numer
            bestEvent = None    #Initialize the event
            minDistance = 1000  #Used to track the minimum distance from singleton to operon that has an identical gene
            
            filteredList = iter(filter(lambda x:x.fragmentIndex == i, fragments)) #Get the fragment we need based on the index
            unmarkedFragment = next(filteredList, None)
            
            if len(unmarkedFragment.sequence) > 1: #We're processing an operon
                for j in range(0, len(coverageTracker)):
                    filteredList = iter(filter(lambda x:x.fragmentIndex == j, fragments)) #Get the fragment we need based on the index
                    currFragment = next(filteredList, None)
                    
                    if i != j and coverageTracker[j] == True and len(currFragment.sequence) > 1:
                        op1 = unmarkedFragment.sequence
                        op2 = currFragment.sequence
                        
                        event = Event(0)
                        event.setFragmentDetails1(unmarkedFragment)
                        event.setFragmentDetails2(currFragment)
                        event.setGenome1Name(strain.name)
                        event.setGenome2Name(strain.name)
                        event.setTechnique('Self Global Alignment (Operon)')
                        event.setAncestralOperonGeneSequence(copy.deepcopy(op1)) #Set the ancestral operon sequence
                        score, event = performGlobalAlignment(op1, op2, event) #Perform the global alignment
                        event.setScore(score)
                        
                        #Compute whether this comparison is interesting
                        threshold = max(len(op1), len(op2))
                        threshold = threshold//3
                        
                        numOperonDifferences = computeOperonDifferences(op1, op2)
                        
                        if numOperonDifferences <= threshold and score < bestScore:
                            bestScore = score
                            bestEvent = event
            else: #We're processing a singleton
                for j in range(0, len(coverageTracker)):
                    if i != j and coverageTracker[j] == True:
                        filteredList = iter(filter(lambda x:x.fragmentIndex == j, fragments)) #Get the fragment we need based on the index
                        currFragment = next(filteredList, None)
                        
                        op1 = unmarkedFragment.sequence
                        op2 = currFragment.sequence
                        
                        event = Event(0)
                        event.setFragmentDetails1(unmarkedFragment)
                        event.setFragmentDetails2(currFragment)
                        event.setGenome1Name(strain.name)
                        event.setGenome2Name(strain.name)
                        event.setTechnique('Self Global Alignment (Singleton)')
                        event.setAncestralOperonGeneSequence(copy.deepcopy(op1)) #Set the ancestral operon sequence
                        
                        if op1[0] in op2 and abs(i-j) < minDistance: #Checks if the singleton gene is located in the operon and if the distance is smaller
                            minDistance = abs(i-j)
                            event.setScore(0)
                            bestEvent = event
                            
            if bestEvent != None: #A match was found meaning the operon is a duplicate therefor do not add it into the ancestor
                globals.trackingId += 1
                bestEvent.trackingEventId = globals.trackingId
                coverageTracker[i] = True
                duplicationEvents.append(bestEvent)
                
                duplicationDetails = ''
                position = bestEvent.fragmentDetails1.startPositionInGenome
                op = copy.deepcopy(bestEvent.fragmentDetails1.sequence)
                
                if bestEvent.fragmentDetails1.isNegativeOrientation == True: #Reverses the genes if the operon was originally negative to ensure the correct position is computed
                    op.reverse()
                
                for gene in op:
                    duplicationDetails += gene + ' ' + str(position) + ','
                    position += 1
                duplicationDetails += '|'
                
                #Increment the duplicate counter with size of operon since the operon is a duplication
                strain = addDuplicationEventsToStrain(strain, [len(bestEvent.fragmentDetails1.sequence)], duplicationDetails)
                
                print('\n&&&&&& Self Global Alignment &&&&&')
                print(bestEvent.toString())
                print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
        
            else: #No match was found therefore it must have been lost in the sibling therefore we include it in the ancestor
                coverageTracker[i] = True
                globals.trackingId +=1
                    
                event = Event(globals.trackingId)    
                event.setScore(-1)
                event.setFragmentDetails1(unmarkedFragment)
                event.setFragmentDetails2(unmarkedFragment)
                event.setGenome1Name(strain.name)
                event.setGenome2Name(strain.name)
                event.setTechnique('Self Global Alignment (No match!)')
                event.setAncestralOperonGeneSequence(copy.deepcopy(unmarkedFragment.sequence))
                lossEvents.append(event)
                
                deletionDetails = ''
                position = event.fragmentDetails1.startPositionInGenome
                op = copy.deepcopy(event.fragmentDetails1.sequence)
                
                if event.fragmentDetails1.isNegativeOrientation == True: #Reverses the genes if the operon was originally negative to ensure the correct position is computed
                    op.reverse()
                
                for gene in op:
                    deletionDetails += gene + ' ' + str(position) + ','
                    position += 1
                deletionDetails += '|'
                
                #Increment the loss counter with the size of the operon since the operon is a loss
                strain = addDeletionEventsToStrain(strain, [len(event.fragmentDetails1.sequence)], deletionDetails)
                
                print('\n&&&&&& Self Global Alignment &&&&&')
                print(event.toString())
                print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
                
    return duplicationEvents, lossEvents, coverageTracker, strain