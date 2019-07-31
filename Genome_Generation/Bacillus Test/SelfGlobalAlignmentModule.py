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
            bestJ = -999
            
            filteredList = iter(filter(lambda x:x.fragmentIndex == i, fragments)) #Get the fragment we need based on the index
            unmarkedFragment = next(filteredList, None)
            
            if len(unmarkedFragment.sequence) > 1: #We're processing an operon
                for j in range(0, len(coverageTracker)):
                    filteredList = iter(filter(lambda x:x.fragmentIndex == j, fragments)) #Get the fragment we need based on the index
                    currFragment = next(filteredList, None)
                    
                    if i != j and currFragment.isDuplicate == False and len(currFragment.sequence) > 1:
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
                        if (score > 0 and score > bestScore) or (score == bestScore and abs(i-j) < minDistance):
                            bestScore = score
                            bestEvent = event
                            minDistance = abs(i-j)
                            bestJ = j
            #Make sure an origin or a terminus doesn't get mapped with a singleton gene
            elif len(unmarkedFragment.sequence) == 1 and unmarkedFragment.description != 'Origin' and unmarkedFragment.description != 'Terminus':
                for j in range(0, len(coverageTracker)):
                    
                    filteredList = iter(filter(lambda x:x.fragmentIndex == j, fragments)) #Get the fragment we need based on the index
                    currFragment = next(filteredList, None)
                    
                    if i != j and currFragment.isDuplicate == False:
                        
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
                            bestJ = j
                            
            if bestEvent != None: #A match was found meaning the operon is a duplicate therefor do not add it into the ancestor
                #Handle special case where two unmarked operons are selected as the best matches
                cycleDuplication = False
                if coverageTracker[i] == False and coverageTracker[bestJ] == False:
                    bestEvent.fragmentDetails2.isDuplicate = True
                    coverageTracker[bestJ] = True
                    cycleDuplication = True
                    bestEvent.setScore(-1)
                    bestEvent.setAncestralOperonGeneSequence(copy.deepcopy(bestEvent.fragmentDetails2.sequence)) #insert source as ancestral operon
                    bestEvent.setTechnique('Self Global Alignment (Cyclic Duplication!)')
                    lossEvents.append(bestEvent)
                    
                globals.trackingId += 1
                coverageTracker[i] = True
                bestEvent.trackingEventId = globals.trackingId
                
                if len(bestEvent.fragmentDetails1.sequence) > 1 and len(bestEvent.fragmentDetails2.sequence) > 1:    
                    handleDuplicateDetails(bestEvent, strain, sibling, cycleDuplication)
                else:
                    #Singleton was mapped to an operon
                    if len(bestEvent.fragmentDetails1.sequence) == 1:
                        gene = bestEvent.fragmentDetails1.sequence[0]
                        position = bestEvent.fragmentDetails1.startPositionInGenome
                    else:
                        gene = bestEvent.fragmentDetails2.sequence[0]
                        position = bestEvent.fragmentDetails2.startPositionInGenome
                    tempString = gene + ' ' + str(position) + ';'
                    strain = addDuplicationEventsToStrain(strain, [1], tempString)
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
                
                position = event.fragmentDetails1.startPositionInGenome
                op = copy.deepcopy(event.fragmentDetails1.sequence)
                
                tempString = ''
                for n in range(0, len(op)):
                    gene = op[n]
                    if event.fragmentDetails1.isNegativeOrientation == False:
                        tempString += gene + ' ' + str(n + position) + ', '
                    else:
                        tempString = gene + ' ' + str(position + len(op) - n - 1) + ', ' + tempString
                tempString = tempString[0:(len(tempString) - 2)] #Remove the last comma and space
                tempString += ';'
                
                #Increment the loss counter with the size of the operon since the operon is a loss
                sibling = addDeletionEventsToStrain(sibling, [len(event.fragmentDetails1.sequence)], tempString)
                
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
def handleDuplicateDetails(event, strain, sibling, cycleDuplication):
    operon1Gaps = event.operon1Gaps
    operon2Gaps = event.operon2Gaps
    operon1GapPositions = event.operon1GapPositions    
    
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
            for m in range(0, len(operon1Gaps)):
                gap = operon1Gaps[m]
                positions = operon1GapPositions[m]

                for n in range(0, len(gap)):
                    gene = gap[n]
                    position = positions[n]
                    if event.fragmentDetails1.isNegativeOrientation: #This compute the correct gene position based on whether operon was in the negative orientation or not originally
                        genePos = event.fragmentDetails1.startPositionInGenome + len(event.fragmentDetails1.sequence) - position - 1
                        tempString = gene + ' ' + str(genePos) + ', ' + tempString
                    else:
                        genePos = position + event.fragmentDetails1.startPositionInGenome
                        tempString += gene + ' ' +str(genePos) + ', '
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
    sequenceDuplicated = event.fragmentDetails2.sequence
    strain.duplicationDetails += event.selfDuplication #Built from the trace back
    
    sizeOfDuplication = len(sequenceDuplicated)
    if sizeOfDuplication in strain.duplicationCounts:
        strain.duplicationCounts[sizeOfDuplication] += 1
    else:
        strain.duplicationCounts[sizeOfDuplication] = 1
        
    #Special case when two unmarked operons are mapped as good matches then we must indicate the operon was lost in the sibling
    if cycleDuplication:
        sibling.deletionDetails += event.selfDuplication
        if sizeOfDuplication in sibling.deletionCounts:
            sibling.deletionCounts[sizeOfDuplication] += 1
        else:
            sibling.deletionCounts[sizeOfDuplication] = 1
            
    ####################################
    #Handle the Codon Mismatches here##
    ###################################
    if len(event.codonMismatchGenesStrain1) > 0 and globals.enableDeletionReversions:
        for w in range(0, len(event.codonMismatchGenesStrain1)):
            temp = event.codonMismatchGenesStrain1[w].split('-')
            gene = temp[0].strip()            
            if event.fragmentDetails1.isNegativeOrientation == False:
                position = event.codonMismatchIndexesStrain1[w] + event.fragmentDetails1.startPositionInGenome
            else:
                position = event.fragmentDetails1.startPositionInGenome + len(event.fragmentDetails1.sequence) - event.codonMismatchIndexesStrain1[w] - 1
            globals.codonMismatchCounter += 1 #Increment the counter by one for each codon mismatch
            strain.tempCodonDetails += gene + ' ' + str(position) + ';'
            
    ################################
    #Handle the Substitutions here##
    ################################
    if len(event.substitutionGenesStrain1) > 0 and globals.enableDeletionReversions:
        for w in range(0, len(event.substitutionGenesStrain1)):
            temp = event.substitutionGenesStrain1[w].split('-')
            gene = temp[0].strip()            
            if event.fragmentDetails1.isNegativeOrientation == False:
                position = event.substitutionIndexesStrain1[w] + event.fragmentDetails1.startPositionInGenome
            else:
                position = event.fragmentDetails1.startPositionInGenome + len(event.fragmentDetails1.sequence) - event.substitutionIndexesStrain1[w] - 1
            globals.substitutionCounter += 1 #Increment the counter by one for each substitution
            strain.tempSubstitutionDetails += gene + ' ' + str(position) + ';'
            
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