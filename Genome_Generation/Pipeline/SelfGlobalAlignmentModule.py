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
def findOrthologsBySelfGlobalAlignment(strain, coverageTracker, sibling):
    if globals.printToConsole:
        print('Performing self-global alignment on strain: %s' %(strain.name))
        
    lossEvents = []
    duplicationEvents = []
    fragments = strain.genomeFragments
    
    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False: #Operon has not been marked
            bestScore = -1000    #Make the best score some large numer
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
                        #threshold = max(len(op1), len(op2))
                        #threshold = threshold//3
                        #numOperonDifferences = computeOperonDifferences(op1, op2)
                        #if numOperonDifferences <= threshold and score < bestScore:
                        if score > 0 and score > bestScore and abs(i-j) < minDistance:
                            bestScore = score
                            bestEvent = event
                            minDistance = abs(i-j)
            #Make sure an origin or a terminus doesn't get mapped with a singleton gene
            elif len(unmarkedFragment.sequence) == 1 and unmarkedFragment.description != 'Origin' and unmarkedFragment.description != 'Terminus':
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
                
                handleDuplicateDetails(bestEvent, strain)
                
                coverageTracker[i] = True
                duplicationEvents.append(bestEvent)
                
                #This is now being handled with the function handleDuplicateDetails
                #duplicationDetails = ''
                #position = bestEvent.fragmentDetails1.startPositionInGenome
                #op = copy.deepcopy(bestEvent.fragmentDetails1.sequence)
                
                #if bestEvent.fragmentDetails1.isNegativeOrientation == True: #Reverses the genes if the operon was originally negative to ensure the correct position is computed
                    #op.reverse()
                
                #for gene in op:
                    #duplicationDetails += gene + ' ' + str(position) + ', '
                    #position += 1
                #duplicationDetails = duplicationDetails[0:(len(duplicationDetails) - 2)]
                #duplicationDetails += ';'
                
                #Increment the duplicate counter with size of operon since the operon is a duplication
                #strain = addDuplicationEventsToStrain(strain, [len(bestEvent.fragmentDetails1.sequence)], duplicationDetails)
                if globals.printToConsole:
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
                
#                deletionDetails = ''
#                position = event.fragmentDetails1.startPositionInGenome
#                op = copy.deepcopy(event.fragmentDetails1.sequence)
#                
#                if event.fragmentDetails1.isNegativeOrientation == True: #Reverses the genes if the operon was originally negative to ensure the correct position is computed
#                    op.reverse()
#                
#                for gene in op:
#                    deletionDetails += gene + ' ' + str(position) + ', '
#                    position += 1
#                deletionDetails = deletionDetails[0:(len(deletionDetails) - 2)]
#                deletionDetails += ';'
#                
#                #Increment the loss counter with the size of the operon since the operon is a loss
#                sibling = addDeletionEventsToStrain(sibling, [len(event.fragmentDetails1.sequence)], deletionDetails)
                if globals.printToConsole:
                    print('\n&&&&&& Self Global Alignment &&&&&')
                    print(event.toString())
                    print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
                
    return duplicationEvents, lossEvents, coverageTracker, strain, sibling

######################################################
# handleDuplicateDetails
# Parameters:
# Description: Reconstructs the ancestral operon by determining whether the gaps are losses or duplications
######################################################
def handleDuplicateDetails(event, strain):
    operon1Gaps = event.operon1Gaps
    operon2Gaps = event.operon2Gaps
    
    if len(operon1Gaps) > 0 or len(operon2Gaps) > 0:
        #S1 is the duplicate operon ie target therefore S2 is the source
        if len(operon2Gaps) > 0:
            #These are the extra genes in the source therefore they were deleted in the target
            tempString = ''
            for gap in operon2Gaps:
                tempString = ''
                for gene in gap:
                    if event.fragmentDetails2.isNegativeOrientation:
                        tempString = '!' + gene + ' ' +str(-1) + ', ' + tempString
                    else:
                        tempString += '!' + gene + ' ' +str(-1) + ', '
                tempString = tempString[0:(len(tempString) - 2)] #Remove the last comma and space
                tempString += ';'
                strain.deletionDetails += tempString
                
                #Increment the size distribution for this particular gap
                sizeOfDeltion = len(gap)
                if sizeOfDeltion in strain.deletionCounts:
                    strain.deletionCounts[sizeOfDeltion] += 1
                else:
                    strain.deletionCounts[sizeOfDeltion] = 1
                    
        if len(operon1Gaps) > 0:
            #These are the extra genes in the target, therefore there were duplicated into the target
            tempString = ''
            for gap in operon1Gaps:
                tempString = ''
                for gene in gap:
                    if event.fragmentDetails1.isNegativeOrientation:
                        tempString = '!' + gene + ' ' +str(-1) + ', ' + tempString
                    else:
                        tempString += '!' + gene + ' ' +str(-1) + ', '
                tempString = tempString[0:(len(tempString) - 2)] #Remove the last comma and space
                tempString += ';'
                strain.duplicationDetails += tempString
                #Increment the size distribution
                sizeOfDuplication = len(gap)
                if sizeOfDuplication in strain.duplicationCounts:
                    strain.duplicationCounts[sizeOfDuplication] += 1
                else:
                    strain.duplicationCounts[sizeOfDuplication] = 1
                    
    #Indicate the whole operon was duplicated
    index = 0
    tempString = ''
    sequenceTarget = event.fragmentDetails1.sequence
    sequenceDuplicated = event.fragmentDetails2.sequence
    position = event.fragmentDetails1.startPositionInGenome
    rememberPoint = 0
    for x in range(0, len(sequenceDuplicated)):
        found = False
        while index < len(sequenceTarget):
            if sequenceTarget[index] == sequenceDuplicated[x]:
                if event.fragmentDetails1.isNegativeOrientation == False:
                    tempString += sequenceDuplicated[x] + ' ' + str(index + position) + ', '
                else:
                    tempString = sequenceDuplicated[x] + ' ' + str(position + len(event.fragmentDetails1.sequence) - index - 1) + ', ' + tempString
                rememberPoint = index + 1
                found = True
                break
            else:
                index+=1
        if found == False:
            if event.fragmentDetails1.isNegativeOrientation == False:
                tempString += '!' + sequenceDuplicated[x] + ' ' + str(-1) + ', '
            else:
                tempString = '!' + sequenceDuplicated[x] + ' ' + str(-1) + ', ' + tempString
            index = rememberPoint #reset the index back to our last found point
        else:
            index = rememberPoint
    tempString = tempString[0:(len(tempString) - 2)] #Remove the last comma and space
    tempString += ';'
    strain.duplicationDetails += tempString
    
    sizeOfDuplication = len(sequenceDuplicated)
    if sizeOfDuplication in strain.duplicationCounts:
        strain.duplicationCounts[sizeOfDuplication] += 1
    else:
        strain.duplicationCounts[sizeOfDuplication] = 1
        
    #Now handle the match region details, the matched regions are duplications with a normal index
#    tempString = ''
#    sequence = event.ancestralOperonGeneSequence
#    position = event.fragmentDetails1.startPositionInGenome
#    for x in range(0, len(sequence)):
#        if event.fragmentDetails1.isNegativeOrientation == False:
#            tempString += sequence[x] + ' ' + str(x + position) + ', '
#        else:
#            tempString = sequence[x] + ' ' + str(position + len(sequence) - x - 1) + ', ' + tempString
#    tempString = tempString[0:(len(tempString) - 2)] #Remove the last comma and space
#    tempString += ';'
#    strain.duplicationDetails += tempString
#    
#    sizeOfDuplication = len(sequence)
#    if sizeOfDuplication in strain.duplicationCounts:
#        strain.duplicationCounts[sizeOfDuplication] += 1
#    else:
#        strain.duplicationCounts[sizeOfDuplication] = 1