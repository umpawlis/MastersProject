import multiset
import globals

####################################
##Sequence Service Functions########
####################################

######################################################
# computeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: Computes the gene content difference between the operons (multiset) ie whats the difference between them while not removing duplicates
######################################################
def computeOperonDifferences(operon1, operon2):
    set1 = multiset.Multiset();
    set2 = multiset.Multiset();
    
    for op in operon1:
        set1.add(op.split('_')[0].strip())
    for op in operon2:
        set2.add(op.split('_')[0].strip())
        
    set3 = set1.symmetric_difference(set2)
    
    return len(set3)

######################################################
# formatAndComputeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: Creates an array of operons 1 and 2, computes how many genes are unique to 1 operon (no duplicates) and computes the gene content difference between the operons (multiset) ie whats the difference between them while not removing duplicates
######################################################
def formatAndComputeOperonDifferences(operon1, operon2):
    noDuplicatesSet1 = set()
    noDuplicatesSet2 = set()

    set1 = multiset.Multiset();
    set2 = multiset.Multiset();

    operon1 = operon1.replace('-', '')
    operon1 = operon1.replace('[', '')
    operon1 = operon1.replace(']', '')

    operon2 = operon2.replace('-', '')
    operon2 = operon2.replace('[', '')
    operon2 = operon2.replace(']', '')

    operon1List = operon1.split(',')
    operon2List = operon2.split(',')

    noWhiteSpaceOperon1List = []
    noWhiteSpaceOperon2List = []

    for op in operon1List:
        set1.add(op.split('_')[0].strip())
        noDuplicatesSet1.add(op.split('_')[0].strip())
        noWhiteSpaceOperon1List.append(op.strip())

    for op in operon2List:
        set2.add(op.split('_')[0].strip())
        noDuplicatesSet2.add(op.split('_')[0].strip())
        noWhiteSpaceOperon2List.append(op.strip())

    set3 = set1.symmetric_difference(set2)
    noDuplicatesSet3 = noDuplicatesSet1.symmetric_difference(noDuplicatesSet2)

    return len(set3), noWhiteSpaceOperon1List, noWhiteSpaceOperon2List, len(noDuplicatesSet3)

######################################################
# reverseSequence
# Parameters: operon - operon to check
# Description: checks if the operon needs to be reversed
######################################################
def reverseSequence(operon):
    if '-' in operon:
        return True
    return False

######################################################
# formatAllOperons
# Parameters: sequence
# Description: Formats all operons to correct orientation.
######################################################
def formatAllOperons(sequence):
    sequenceList = []
    operonIndexConversions = []
    operonIndex = 0

    for operon1 in sequence:
        reverseOp = reverseSequence(operon1)

        operon1 = operon1.replace('-', '')
        operon1 = operon1.replace('[', '')
        operon1 = operon1.replace(']', '')

        operon1List = operon1.split(',')

        if len(operon1List) > 1:
            noWhiteSpaceOperon1List = []

            for op in operon1List:
                noWhiteSpaceOperon1List.append(op.strip())

            if reverseOp:
                noWhiteSpaceOperon1List.reverse()
            sequenceList.append(noWhiteSpaceOperon1List)
            operonIndexConversions.append(operonIndex)
            operonIndex += 1
        else:
            operonIndexConversions.append(-1)

    return sequenceList, operonIndexConversions

######################################################
# findUniqueGenes
# Parameters: geneList, sequence
# Description: Tries to find the list of genes in another operon of the genome.
######################################################
def findUniqueGenes(geneList, sequence, opIndex, gapPositions, strainDup, strainDel, operonDupPosition, operonDelPosition, alignedRange=(0,0)):
    comparisonSize = len(geneList)
    startIndex = 0
    currentIndex = 0
    genesNotFound = True
    endOfList = False
    missedGene = False
    newSet = False
    geneRanges = []
    numGeneMatches = 0
    duplicationSizes = []
    updatedGeneList = []

    while (comparisonSize >= 2) and genesNotFound:
        #Slide the window across the gene list, one gene at a time
        while (startIndex < comparisonSize) and (startIndex+comparisonSize <= len(geneList)) and genesNotFound:

            currentIndex = startIndex
            #Slide the window across the gene list by comparison size
            while currentIndex < len(geneList):
                if currentIndex+comparisonSize <= len(geneList):
                    gene = geneList[currentIndex:currentIndex+comparisonSize]
                    newSet = True
                # elif not missedGene:
                #     gene = geneList[currentIndex:len(geneList)]
                #     newSet = True

                if newSet:
                    if geneInSequence(gene, sequence, len(gene), opIndex, alignedRange):
                        newRange = (currentIndex, currentIndex+comparisonSize-1)
                        if not checkOverlap(newRange, geneRanges):
                            geneRanges.append(newRange)
                            numGeneMatches += len(gene)
                            duplicationSizes.append(len(gene))
                            
                            #Position of duplication
                            duplicationLocation = gapPositions[newRange[0]:(newRange[1] + 1)]
                            for num in duplicateLocation:
                                strainDup.duplicationInfo+= str(num + operonDupPosition) + ' '
                            strainDup.duplicationInfo+= ';'
                            
                            # updateDuplicationCounter(len(gene))
                        if numGeneMatches == len(geneList):
                            genesNotFound = False
                    else:
                        missedGene = True

                currentIndex += comparisonSize
                newSet = False

            missedGene = False
            startIndex += 1

        currentIndex = 0
        startIndex = 0
        endOfList = False
        #Decrement the window size
        comparisonSize -= 1

    deletionSizes = []
    inSet = True
    deletionSize = 0

    for index in reversed(range(len(geneList))):
        if not checkOverlap((index,index), geneRanges):
            if not inSet:
                inSet = True
            deletionSize += 1
            strainDel.deletionInfo+= str(gapPositions[index] + operonDelPosition) + ' '
            updatedGeneList.insert(0, geneList[index])
        else:
            geneList.pop(index)
            inSet = False
            if deletionSize != 0:
                deletionSizes.append(deletionSize)
                deletionSize = 0
            updatedGeneList.insert(0, '-')
    strainDel.deletionInfo+=';'
    
    #Special case if all genes in list are losses
    if deletionSize != 0:
        deletionSizes.append(deletionSize)
        deletionSize = 0
        
        for num in gapPositions:
            strainDel.deletionInfo+= str(num + operonDelPosition) + ' '
        strainDel.deletionInfo+= ';'
        
    return numGeneMatches, deletionSizes, duplicationSizes, updatedGeneList, strainDup, strainDel

######################################################
# geneInSequence
# Parameters: gene, sequence, comparisonSize, opIndex
# Description: Determines if the list of genes appear somewhere else in the operon.
######################################################
def geneInSequence(gene, sequence, comparisonSize, opIndex, alignedRange):
    geneFound = False
    currentOpIndex = 0
    currentIndex = 0
    searchOperon = True
    checkRange = False
    rangeList = []
    rangeList.append(alignedRange)

    while currentOpIndex < len(sequence) and not geneFound:
        checkRange = False
        operon = sequence[currentOpIndex]
        if currentOpIndex != opIndex:
            searchOperon = True
        else:
            if alignedRange[0] == 0 and alignedRange[1] == 0:
                searchOperon = False
            else:
                searchOperon = True
                checkRange = True

        if searchOperon:
            while currentIndex+comparisonSize <= len(operon) and not geneFound:
                checkGene = True
                if checkRange:
                    if (currentIndex < alignedRange[0]) or ((currentIndex+comparisonSize-1) > alignedRange[1]):
                        checkGene = False

                if checkGene:
                    operonGene = operon[currentIndex:currentIndex+comparisonSize]
                    if operonGene == gene:
                        geneFound = True

                currentIndex += 1
            currentIndex = 0
        currentOpIndex += 1

    return geneFound

######################################################
# checkOverlap
# Parameters: newRange, rangeList
# Description: Checks if the newly matched genes were already matched within another set.
######################################################
def checkOverlap(newRange, rangeList):
    overlap = False

    for ranges in rangeList:
        if (newRange[0] <= ranges[1]) and (ranges[0] <= newRange[1]):
            overlap = True

    return overlap

######################################################
# addDuplicationEventsToStrain
# Parameters: duplicationSizes, strain
# Description: Increments the duplicate counter for the appropriate strain
######################################################
def addDuplicationEventsToStrain(duplicationSizes, strain):
    if duplicationSizes != None and len(duplicationSizes) > 0:
        for x in range(0, len(duplicationSizes)):
            if duplicationSizes[x] in strain.duplicationSizes:
                strain.duplicationSizes[duplicationSizes[x]] += 1
            else:
                strain.duplicationSizes[duplicationSizes[x]] = 1
    return strain
                
######################################################
# addDeletionEventsToStrain
# Parameters: deletionSizes, strain
# Description: Increments the deletion counter for the appropriate strain
######################################################
def addDeletionEventsToStrain(deletionSizes, strain):
    if deletionSizes != None and len(deletionSizes) > 0:
        for x in range(0, len(deletionSizes)):
            if deletionSizes[x] in strain.deletionSizes:
                strain.deletionSizes[deletionSizes[x]] += 1
            else:
                strain.deletionSizes[deletionSizes[x]] = 1
    return strain
                
######################################################
# incrementDuplicateSizeCounters
# Parameters:
# Description: Increments the counters for the duplicate size counters
######################################################
def incrementDuplicateSizeCounters(duplicationSizes):
    if duplicationSizes != None and len(duplicationSizes) > 0:
        for x in range(0, len(duplicationSizes)):
            #Increment the counter for the entire phylogeny
            if duplicationSizes[x] in globals.sizeDuplications:
                globals.sizeDuplications[duplicationSizes[x]] += 1
            else:
                globals.sizeDuplications[duplicationSizes[x]] = 1

            #Increment the local duplicate counter
            if duplicationSizes[x] in globals.localSizeDuplications:
                globals.localSizeDuplications[duplicationSizes[x]] += 1
            else:
                globals.localSizeDuplications[duplicationSizes[x]] = 1

######################################################
# incrementDeletionSizeCounters
# Parameters:
# Description: Increments the counters for the deletion size counters
######################################################
def incrementDeletionSizeCounters(deletionSizes):
    if deletionSizes != None and  len(deletionSizes) > 0:
        for x in range(0, len(deletionSizes)):
            #Increment the deletion counter for the entire phylogeny
            if deletionSizes[x] in globals.sizeDeletions:
                globals.sizeDeletions[deletionSizes[x]] += 1
            else:
                globals.sizeDeletions[deletionSizes[x]] = 1

            #Increment the local duplicate counter
            if deletionSizes[x] in globals.localSizeDeletions:
                globals.localSizeDeletions[deletionSizes[x]] += 1
            else:
                globals.localSizeDeletions[deletionSizes[x]] = 1