from Bio import Phylo
import os
import sys
import argparse
import random
import copy

aminoAcids = ["Ala_GCU", "Ala_GCC", "Ala_GCA", "Ala_GCG", "Arg_CGU", "Arg_CGC", "Arg_CGA", "Arg_CGG", "Arg_AGA", 
    "Arg_AGG", "Asp_GAU", "Asp_GAC", "Cys_UGU", "Cys_UGC", "Glu_GAA", "Glu_GAG", "Gln_CAA", "Gln_CAG", "Gly_GGU", "Gly_GGC",
    "Gly_GGA", "Gly_GGG", "His_CAU", "HIs_CAC", "Ile_AUU", "Ile_AUC", "Ile_AUA", "Leu_UUA", "Leu_UUG", "Leu_CUU", "Leu_CUC", "Leu_CUA", "Leu_CUG",
    "Lys_AAA", "Lys_AAG", "Met_AUG", "Phe_UUU", "Phe_UUC", "Pro_CCU", "Pro_CCC", "Pro_CCA", "Pro_CCG", "Ser_UCU", "Ser_UCC", "Ser_UCA", "Ser_UCG", "Ser_AGU", "Ser_AGC",
    "Thr_ACU", "Thr_ACC", "Thr_ACA", "Thr_ACG", "Trp_UGG", "Tyr_UAU", "Tyr_UAC", "Val_GUU", "Val_GUC","Val_GUA","Val_GUG"]

genomeDirName = "treeGenomes"
beforeTerminus = []
afterTerminus = []
probDup = 0.0
dup_pValue = 0.0
probLoss = 0.0
loss_pValue = 0.0
probInv = 0.0
inv_pValue = 0.0
probSub = 0.0
probTrans = 0.0
trans_pValue = 0.0

ancestorCounter = 1

totalLosses = {}
totalDuplications = {}
distribInversions = {}
distribTranspositions = {}
distribInvertedTrans = {}
totalInversions = 0
totalTranspositions = 0
totalInvertedTrans = 0

testFolder = ""
equalEvents = False

class Event:
    def __init__(self, eventType, indexRange, genes, prevEventRange = None, operonFormat = None):
        self.type = eventType
        self.range = indexRange
        self.genes = genes
        self.prevEventRange = prevEventRange
        self.operonFormat = operonFormat

    def __repr__(self):
        return "Event(Type = %r, Range = %r, Genes: %r" % (self.type, self.range, self.genes)

    def __str__(self):
        event = ""

        if self.type == "I":
            event = self.operonFormat + "|"
        else:
            # for gene, index in zip(self.genes, self.range):
            #     event += gene + " " + str(index) + ","
            event = ', '.join(gene + " " + str(index) for gene, index in zip(self.genes, self.range))
            event += ";"
            
            if self.type == "T" or self.type == "IT":
                event += "|"

        return event

class Node:
    def __init__(self, before, after, parent, dupEvents, lossEvents, invEvents, subEvents, transEvents, invTransEvents, branchEvents, lineageEvents, inversionBefores):
        self.name = None
        self.before = before
        self.after = after
        self.parent = parent
        self.children = []
        self.dupEvents = dupEvents
        self.lossEvents = lossEvents
        self.invEvents = invEvents
        self.subEvents = subEvents
        self.transEvents = transEvents
        self.invTransEvents = invTransEvents
        self.branchEvents = branchEvents
        self.lineageEvents = lineageEvents
        self.inversionBefores = inversionBefores

    def getGenome(self):
        return self.genome

    def getParent(self):
        return self.parent

    def getChildren(self):
        return self.children

    def getBranchEvents(self):
        return self.branchEvents

    def printNode(self):
        outputFile = open(testFolder + "/generatorOutput.txt", "a+")
        outputFile.write("Strain:%s\n" % (self.name))
        outputFile.write("Codon Mismatch:\n") # % (self.codonMismatch)
        outputFile.write("Substitution:%s\n" % (''.join(str(e) for e in self.subEvents)))
        outputFile.write("Duplication:%s\n" % (''.join(str(e) for e in self.dupEvents)))
        outputFile.write("Deletion:%s\n" % (''.join(str(e) for e in self.lossEvents)))
        outputFile.write("Inversion:%s\n" % (''.join(str(e) for e in self.invEvents)))
        outputFile.write("Transposition:%s\n" % (''.join(str(e) for e in self.transEvents)))
        outputFile.write("Inverted Transposition:%s\n" % (''.join(str(e) for e in self.invTransEvents)))

        # for k, v in sorted(totalLosses.items()):
        #     outputFile.write(" size: %s count: %s," % (str(k), str(v)))
        # for event in self.branchEvents:
        #     print event
        outputFile.close()
# Traverse the Newick tree and add events everytime you find more clades on a branch.
def generateTests(testSetDir, treeStructure, max_length, num_operons, num_events, dupProb, dup_p, lossProb, loss_p, invProb, inv_p, subProb, transProb, trans_p, equal):
    global probDup
    global dup_pValue
    global probLoss
    global loss_pValue
    global probInv
    global inv_pValue
    global probSub
    global probTrans
    global trans_pValue
    global ancestorCounter
    
    global testFolder
    testFolder = testSetDir
    global equalEvents
    equalEvents = equal

    probDup = dupProb
    dup_pValue = dup_p
    probLoss = lossProb
    loss_pValue = loss_p
    probInv = invProb
    inv_pValue = inv_p
    probSub = subProb
    probTrans = transProb
    trans_pValue = trans_p

    if not os.path.exists(testSetDir):
        os.makedirs(testFolder)

    createAncestor(max(max_length, 15), num_operons)
#    print "Ancestor:"
#    print formatGenome(beforeTerminus, afterTerminus)
#    print ""
    sequenceFile = open(testFolder + "/root.txt", "w+")
    sequenceFile.write(formatGenome(beforeTerminus, afterTerminus))
    sequenceFile.close()
    resetCounters()

    root = Node(beforeTerminus, afterTerminus, None, "", "", "", "", "", "", None, None, None)

#    print('Reading in newick tree structure from file: %s...' % (treeStructure))
    newickTree = Phylo.read(treeStructure, 'newick')
    currNode = newickTree.clade

    createFile(testFolder + '/' + "generatorOutput.txt", newickTree)

    listOfEvents = []
    if len(currNode.clades) > 0:
        left = buildTreeData(currNode.clades[0], beforeTerminus, afterTerminus, num_events, listOfEvents, root)
        hasRight = False
        if len(currNode.clades) > 1:
            if len(left.invEvents) > 0:
                # Used so that there is a lower chance of both children having inversions (not common)
                invMultiplier = 0.5
            else:
                invMultiplier = 1.0
            right = buildTreeData(currNode.clades[1], beforeTerminus, afterTerminus, num_events, listOfEvents, root, invMultiplier)
            hasRight = True

            leftEvents = copy.deepcopy(left.branchEvents)
            adjustLossIndexes(left.lossEvents, right.branchEvents, right.inversionBefores, None)
            adjustLossIndexes(right.lossEvents, leftEvents, left.inversionBefores, None)

        printTree(left)
        if hasRight:
            printTree(right)
            
    printTotals()
    ancestorCounter = 1
    
def resetCounters():
    global totalLosses
    global totalDuplications
    global distribInversions
    global distribTranspositions
    global distribInvertedTrans
    global totalInversions
    global totalTranspositions
    global totalInvertedTrans
    
    totalLosses = {}
    totalDuplications = {}
    distribInversions = {}
    distribTranspositions = {}
    distribInvertedTrans = {}
    totalInversions = 0
    totalTranspositions = 0
    totalInvertedTrans = 0

def printTree(currNode):
    if len(currNode.children) > 0:
        printTree(currNode.children[0])
        if len(currNode.children) > 1:
            printTree(currNode.children[1])
    currNode.printNode()


def printTotals():
    outputFile = open(testFolder + "/generatorOutput.txt", "a+")
    outputFile.write("Total Deletions:%s\n" % (','.join(" size: %s count: %s" % (str(k), str(v)) for k, v in sorted(totalLosses.items()))))
    outputFile.write("Total Duplications:%s\n" % (','.join(" size: %s count: %s" % (str(k), str(v)) for k, v in sorted(totalDuplications.items()))))
    outputFile.write("Size Distribution of Inversions:%s\n" % (','.join(" size: %s count: %s" % (str(k), str(v)) for k, v in sorted(distribInversions.items()))))
    outputFile.write("Size Distribution of Transpositions:%s\n" % (','.join(" size: %s count: %s" % (str(k), str(v)) for k, v in sorted(distribTranspositions.items()))))
    outputFile.write("Size Distribution of Inverted Transpositions:%s\n" % (','.join(" size: %s count: %s" % (str(k), str(v)) for k, v in sorted(distribInvertedTrans.items()))))

    outputFile.write("Total Inversions: %s\n" % (str(totalInversions)))
    outputFile.write("Total Transpositions: %s\n" % (str(totalTranspositions)))
    outputFile.write("Total Inverted Transpositions: %s\n" % (str(totalInvertedTrans)))
    outputFile.close()

def buildTreeData(node, before, after, numEvents, events, parent, invMultiplier = 1.0):
    global ancestorCounter
    global totalLosses
    global totalDuplications
    global distribInversions
    global distribTranspositions
    global distribInvertedTrans
    global totalInversions
    global totalTranspositions
    global totalInvertedTrans

    currentBefore = copy.deepcopy(before)
    currentAfter = copy.deepcopy(after)
    currentEvents = copy.deepcopy(events)
    branchEvents = []

    duplicationEvents = []
    lossEvents = []
    inversionEvents = []
    substitutionEvents = []
    transpositionEvents = []
    invertedTransEvents = []
    prevEventRange = []
    sectionReversed = False
    inversionBefores = []

    originalIndexes = calculateIndexes(currentBefore, currentAfter)
    currentIndexes = copy.deepcopy(originalIndexes)

#    print "Previous Events:"
#    print currentEvents
#    print "Parent Genome:"
#    print formatGenome(before, after)
#    print ""

    count = 0
    rand = 0.10
    while count < numEvents:
        newEvent = True
        
        if not equalEvents:
            rand = random.random()
        else:
            rand = (rand + 0.25) % 1.0

        if rand < probDup:
            event, eventRange, genes, prevEventRange, sectionReversed = performDuplication(currentBefore, currentAfter, dup_pValue)
            duplicationEvents.append(event)

            if len(genes) in totalDuplications:
                totalDuplications[len(genes)] += 1
            else:
                totalDuplications[len(genes)] = 1
        elif rand < probDup + probLoss:
            event, eventRange, genes = performLoss(currentBefore, currentAfter, loss_pValue)
            lossEvents.append(event)

            if len(genes) in totalLosses:
                totalLosses[len(genes)] += 1
            else:
                totalLosses[len(genes)] = 1
        elif rand < probDup + probLoss + probInv:
            if not equalEvents:
                if random.random() < invMultiplier:
                    event, eventRange, genes, size = performInversion(currentBefore, currentAfter, inv_pValue)
                    inversionEvents.append(event)
                    invMultiplier *= 0.5
                else:
                    newEvent = False
                    count -= 1
            else:
                event, eventRange, genes, size = performInversion(currentBefore, currentAfter, inv_pValue)
                inversionEvents.append(event)
                
            if newEvent:
                if size in distribInversions:
                    distribInversions[size] += 1
                else:
                    distribInversions[size] = 1
                totalInversions += 1
                
        elif rand < probDup + probLoss + probInv + probSub:
            event, eventRange, genes = performSubstitution(currentBefore, currentAfter)
            substitutionEvents.append(event)

        elif rand < probDup + probLoss + probInv + probSub + probTrans:
            event, eventRange, genes, prevEventRange, sectionReversed = performTransposition(currentBefore, currentAfter, trans_pValue)
            if event.type == "T":
                transpositionEvents.append(event)
    
                if len(genes) in distribTranspositions:
                    distribTranspositions[len(genes)] += 1
                else:
                    distribTranspositions[len(genes)] = 1
                totalTranspositions += 1
            elif event.type == "IT":
                invertedTransEvents.append(event)
                
                if len(genes) in distribInvertedTrans:
                    distribInvertedTrans[len(genes)] += 1
                else:
                    distribInvertedTrans[len(genes)] = 1
                totalInvertedTrans += 1

        if newEvent:
            currentEvents.append(event)
            if event.type == "I":
                inversionBefores.append(copy.deepcopy(currentBefore))

#            print branchEvents
#            print event
            updateEvents(branchEvents, event, currentBefore, currentAfter, prevEventRange)
            updateIndexes(event, originalIndexes, currentIndexes, prevEventRange, sectionReversed)
            branchEvents.append(event)
#            print branchEvents

        count += 1

    currNode = Node(currentBefore, currentAfter, parent, duplicationEvents, lossEvents, inversionEvents, substitutionEvents, transpositionEvents, invertedTransEvents, branchEvents, currentEvents, inversionBefores)

    if len(node.clades) > 0:
        left = buildTreeData(node.clades[0], currentBefore, currentAfter, numEvents, currentEvents, currNode)
        currNode.children.append(left)
        hasRight = False
        if len(node.clades) > 1:
            if len(left.invEvents) > 0:
                # Used so that there is a lower chance of both children having inversions (not common)
                invMultiplier = 0.5
            else:
                invMultiplier = 1.0
            right = buildTreeData(node.clades[1], currentBefore, currentAfter, numEvents, currentEvents, currNode, invMultiplier)
            currNode.children.append(right)
            hasRight = True

            # print "Updating Loss Events"
            # print left.lossEvents
            # print right.lossEvents

            leftEvents = copy.deepcopy(left.branchEvents)
            adjustLossIndexes(left.lossEvents, right.branchEvents, right.inversionBefores, after)
            adjustLossIndexes(right.lossEvents, leftEvents, left.inversionBefores, after)

            # print left.lossEvents
            # print right.lossEvents

        currNode.name = "Ancestor " + str(ancestorCounter)
        sequenceFile = open(testFolder + "/" + currNode.name + ".txt", "w+")
        sequenceFile.write(formatGenome(currentBefore, currentAfter))
        sequenceFile.close()
        ancestorCounter += 1
        # currNode.printNode()
    elif node.name is not None and len(node.name) > 0:
        completePath = testFolder + "/" + node.name
        if not os.path.exists(completePath):
            os.makedirs(completePath)

        sequenceFile = open(completePath + "/sequence.txt", "w+")
        sequenceFile.write(formatGenome(currentBefore, currentAfter))
        sequenceFile.close()
        eventsFile = open(completePath + "/events.txt", "w+")
        for eachEvent in currentEvents:
            eventsFile.write(repr(eachEvent) + "\n")
        eventsFile.close()

        currNode.name = node.name
        # currNode.printNode()

    return currNode

def calculateIndexes(before, after):
    count = 2

    for item in before:
        if isinstance(item, list):
            count += len(item)
        else:
            count += 1
    for item in after:
        if isinstance(item, list):
            count += len(item)
        else:
            count += 1

    return range(count)

def updateIndexes(event, origIndexes, currIndexes, prevEventRange = None, sectionReversed = None):
    overlap, overlapRange = checkOverlap(currIndexes, event.range, event.type)

    if event.type == "D":
        for index1 in prevEventRange:
            origIndexes.append(origIndexes[currIndexes.index(index1)])
        for i in range(len(currIndexes)):
            if currIndexes[i] >= event.range[0]:
                currIndexes[i] = currIndexes[i] + len(event.range)

        copyRange = range(event.range[0], event.range[-1]+1)
        if sectionReversed:
            copyRange.reverse()

        for index2 in copyRange:
            currIndexes.append(index2)
    elif event.type == "I":
        ordered = getOrderedIndexes(currIndexes, event.range)
        copyRange = range(event.range[0], event.range[-1]+1)
        copyRange.reverse()

        for index1, index2 in zip(ordered, copyRange):
            currIndexes[index1] = index2
    elif event.type == "T" or event.type == "IT":
#        print prevEventRange
        ordered = getOrderedIndexes(currIndexes, prevEventRange)
        for i in range(len(currIndexes)):
            if currIndexes[i] > prevEventRange[-1]:
                currIndexes[i] = currIndexes[i] - len(prevEventRange)
            if currIndexes[i] >= event.range[0]:
                currIndexes[i] = currIndexes[i] + len(event.range)

        if sectionReversed:
            ordered.reverse()

        for index1, index2 in zip(ordered, event.range):
            currIndexes[index1] = index2
    elif event.type == "L":
        ordered = getOrderedIndexes(currIndexes, event.range)
#        print event.range[-1]
        for i in range(len(currIndexes)):
            if currIndexes[i] > event.range[-1]:
                # print "Index " + str(i)
                # print currIndexes[i]
                currIndexes[i] = currIndexes[i] - len(event.range)
                # print currIndexes[i]

        for index in overlap:
            currIndexes[index] = -1

        for i, index in zip(range(len(event.range)), ordered):
            event.range[i] = origIndexes[index]

#    print origIndexes
#    print currIndexes

def getOrderedIndexes(currIndexes, eventRange):
    orderedIndexes = []

    for index in range(eventRange[0], eventRange[-1]+1):
        for index2 in range(len(currIndexes)):
            if index == currIndexes[index2]:
                orderedIndexes.append(index2)
                break

    return orderedIndexes

def updateEvents(eventList, newEvent, before, after, prevEventRange = None, updateLosses = False):
    for event in eventList:
#        print newEvent.range
#        print event
        overlap, overlapRange = checkOverlap(event.range, newEvent.range, newEvent.type)

        if event.type != "L" or updateLosses:
            if newEvent.type == "D":
                # if checkBeforeOrOverlap(event.range, newEvent.range):
                for i in range(len(event.range)):
                    if event.range[i] >= newEvent.range[0]:
                        event.range[i] = event.range[i] + len(newEvent.range)
            elif newEvent.type == "I":
                terminusIndex = getTerminusIndex(before)
                numBefore = 0
                numAfter = 0

                for index in newEvent.range:
                    if index > terminusIndex:
                        numAfter += 1
                    else:
                        numBefore += 1

                oldterminusIndex = terminusIndex + (numAfter - numBefore)
#                print "Terminus index: %d" % (oldterminusIndex)
                for index in overlap:
                    diff = abs(oldterminusIndex - event.range[index])
#                    print "index: %d diff: %d" % (event.range[index], diff)
                    if event.range[index] > oldterminusIndex:
                        event.range[index] = terminusIndex - diff
                    else:
                        event.range[index] = terminusIndex + diff
            elif newEvent.type == "S":
                if overlap and (len(overlap) == 1):
                    event.genes[overlap[0]] = newEvent.genes[0]
            elif newEvent.type == "T" or newEvent.type == "IT":
                prevOverlap, overlapRange = checkOverlap(prevEventRange, event.range, newEvent.type)
                eventOverlap, overlapRange = checkOverlap(event.range, prevEventRange, newEvent.type)
                for i in range(len(event.range)):
                    if event.range[i] > prevEventRange[-1]:
                        event.range[i] = event.range[i] - len(prevEventRange)
                    if event.range[i] >= newEvent.range[0]:
                        if event.range[i] not in overlapRange:
                            event.range[i] = event.range[i] + len(newEvent.range)
                for index in eventOverlap:
                    for index2 in prevOverlap:
#                        print "Prev Event Range"
#                        print prevEventRange
#                        print "Event range"
#                        print event.range
#                        print "Prev Event overlap"
#                        print prevOverlap
#                        print "Event overlap"
#                        print eventOverlap
                        if event.range[index] == prevEventRange[index2]:
                            event.range[index] = newEvent.range[index2]
                            break
            elif newEvent.type == "L":
                for i in range(len(event.range)):
                    if event.range[i] > newEvent.range[-1]:
                        event.range[i] = event.range[i] - len(newEvent.range)
                for index in overlap:
                    event.range[index] = -1


def checkOverlap(range1, range2, eventType):
    overlapIndexes = []
    overlapRange = []

    for i in range(len(range1)):
        if range1[i] in range2:
#            if eventType == "T":
#                print range1[i]
#                print range2
            overlapIndexes.append(i)
            overlapRange.append(range1[i])
        elif eventType == "I":
            if ((range2[0] <= range1[i]) and (range2[-1] >= range1[i])):
                overlapIndexes.append(i)
                overlapRange.append(range1[i])

    return overlapIndexes, overlapRange

def getTerminusIndex(before):
    index = 1

    for item in before:
        if isinstance(item, list):
            index += len(item)
        else:
            index += 1
    return index

def checkBeforeOrOverlap(range1, range2):
    return (((range1[0] <= range2[0]) and (range1[-1] >= range2[0])) or (range1[0] > range2[0])) 

def adjustLossIndexes(losses, branchEvents, befores, after):
    lossUpdate = True
    count = 0

    for event in branchEvents:
        if event.type == "I":
            updateEvents(losses, event, befores[count], after, event.prevEventRange, lossUpdate)
            count += 1
        else:
            updateEvents(losses, event, None, after, event.prevEventRange, lossUpdate)

def formatGenome(before, after):
    sequence = "< o >, "

    for i in range(0, len(before)):
        if isinstance(before[i], list):
            sequence += "["
            for j in range(0, len(before[i])):
                sequence += before[i][j]
                if j != len(before[i])-1:
                    sequence += ", "
            sequence += "]"
        else:
            sequence += before[i]

        if i != len(before)-1:
            sequence += ", "

    sequence += ", < t >, "

    for i in range(0, len(after)):
        if isinstance(after[i], list):
            sequence += "-["
            for j in range(0, len(after[i])):
                sequence += after[i][j]
                if j != len(after[i])-1:
                    sequence += ", "
            sequence += "]"
        else:
            sequence += "-" + after[i]

        if i != len(after)-1:
            sequence += ", "

    return sequence


def performTransposition(before, after, p):
    event = ""
    genes = []
    startPos = 0
    operonTrans = False
    sectionReversed = False

    if random.random() < 0.5:
        genome = before
        fromBefore = True
    else:
        genome = after
        fromBefore = False

    index = random.randint(0, len(genome)-1)
    if isinstance(genome[index], list):
        # print "is list"
        operon = genome[index]
        # stopPos = min(geometricSampling(p), len(operon))
        stopPos = len(operon)
        transSection = operon[startPos:stopPos]

        if len(transSection) == len(operon):
            del genome[index]
            operonTrans = True
        else:
            del genome[index][startPos:stopPos]

    else:
        # print "is singleton"
        transSection = genome[index:index+1]
        del genome[index]

    beforeStartIndex = getAbsoluteIndex(fromBefore, before, after, index, startPos)
    beforeStopIndex = beforeStartIndex + len(transSection)

    if random.random() < 0.5:
        genome = before
        targetBefore = True
    else:
        genome = after
        targetBefore = False
        
#    print transSection
    eventType = "T"
    if fromBefore and not targetBefore:
        transSection.reverse()
        sectionReversed = True
        eventType = "IT"
    elif not fromBefore and targetBefore:
        transSection.reverse()
        sectionReversed = True
        eventType = "IT"

    targetIndex = random.randint(0, len(genome)-1)
    targetPos = 0
    
    if isinstance(genome[targetIndex], list):
        # print "target in list"
        operon = genome[targetIndex]
        targetPos = random.randint(0, len(operon)-1)

        for gene in reversed(transSection):
            # print "inserting " + gene
            genome[targetIndex].insert(targetPos, gene)

        # event = "Performing transposition... From before: %s Index: %s  Section length: %s Target before: %s Target index: %s Target position: %s" % (str(fromBefore), str(index), str(len(transSection)), str(targetBefore), str(targetIndex), str(targetPos))
#        print "Performing transposition... From before: %s Index: %s  Section length: %s Target before: %s Target index: %s Target position: %s" % (str(fromBefore), str(index), str(len(transSection)), str(targetBefore), str(targetIndex), str(targetPos))
    else:
        # print "target in singleton"
        if operonTrans:
            genome.insert(targetIndex, transSection)
        else:
            for gene in reversed(transSection):
                # print "inserting " + gene
                genome.insert(targetIndex, gene)

        # event = "Performing transposition... From before: %s Index: %s Section length: %s Target before: %s Target index: %s" % (str(fromBefore), str(index), str(len(transSection)), str(targetBefore), str(targetIndex))
#        print "Performing transposition... From before: %s Index: %s Section length: %s Target before: %s Target index: %s" % (str(fromBefore), str(index), str(len(transSection)), str(targetBefore), str(targetIndex))

    absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)
    startIndex = absoluteIndex

    # if absoluteIndex < beforeStartIndex:
    #     beforeStartIndex += len(transSection)
    #     beforeStopIndex += len(transSection)

    for gene in transSection:
        event += gene + " " + str(absoluteIndex) + ","
        genes.append(gene)
        absoluteIndex += 1
#    print event
    stopIndex = absoluteIndex

#    print "Transposition:"
#    print formatGenome(before, after)
#    print ""

    branchEvent = Event(eventType, range(startIndex, stopIndex), genes, range(beforeStartIndex, beforeStopIndex))
    return branchEvent, range(startIndex, stopIndex), genes, range(beforeStartIndex, beforeStopIndex), sectionReversed

def performSubstitution(before, after):
    gene = []
    if random.random() < 0.5:
        genome = before
        fromBefore = True
    else:
        genome = after
        fromBefore = False

    index = random.randint(0, len(genome)-1)
    if isinstance(genome[index], list):
        # print "is list"
        targetPos = random.randint(0, len(genome[index])-1)

        genome[index][targetPos] = random.choice(aminoAcids)
        geneSub = genome[index][targetPos]
    else:
        # print "is singleton"
        targetPos = 0
        genome[index] = random.choice(aminoAcids)
        geneSub = genome[index]

    absoluteIndex = getAbsoluteIndex(fromBefore, before, after, index, targetPos)
    gene.append(geneSub)

    event = "%s %d" % (geneSub, absoluteIndex)
#    print "Performing substitution... From before: %s Index: %s " % (str(fromBefore), str(index))
#    print "Substitution:"
#    print formatGenome(before, after)
#    print ""

    eventType = "S"
    branchEvent = Event(eventType, range(absoluteIndex, absoluteIndex+1), gene)
    return branchEvent, range(absoluteIndex, absoluteIndex+1), gene

def getDifferentGene(currGene):
    while True:
        newGene = random.choice(aminoAcids)

        if currGene != newGene:
            break

    return newGene

def performInversion(before, after, p):
    event = ""
    genes = []
    lengthBefore = min(geometricSampling(p, 0), len(before))
    if not lengthBefore:
        lengthAfter = min(geometricSampling(p), len(after))
    else:
        lengthAfter = min(geometricSampling(p, 0), len(after))

    beforeSection = before[len(before)-lengthBefore:]
    afterSection = after[:lengthAfter]
    index = len(before)-lengthBefore
    # print len(before)
    # print before
    # print lengthBefore
    fromBefore = True

    del before[len(before)-lengthBefore:]
    del after[:lengthAfter]
    if afterSection:
        reverseOperons(afterSection)
        before.extend(afterSection[::-1])
    if beforeSection:
        reverseOperons(beforeSection)
        after[0:0] = beforeSection[::-1]

    indexes = []
    absoluteIndex = getAbsoluteIndex(fromBefore, before, after, index, 0)

    if afterSection:
        for operon in reversed(afterSection):
            if isinstance(operon, list):
                for gene in operon:
                    event += gene + " " + str(absoluteIndex) + ", "
                    genes.append(gene)
                    indexes.append(absoluteIndex)
                    absoluteIndex += 1
                event = event[:-2]
                event += ";"
            else:
                event += operon + " " + str(absoluteIndex) + ";"
                genes.append(operon)
                indexes.append(absoluteIndex)
                absoluteIndex += 1

    event += "< t > " + str(absoluteIndex) + ";"
    absoluteIndex += 1 # for the terminus
    if beforeSection:
        for operon in reversed(beforeSection):
            if isinstance(operon, list):
                for gene in operon:
                    event += gene + " " + str(absoluteIndex) + ", "
                    genes.append(gene)
                    indexes.append(absoluteIndex)
                    absoluteIndex += 1
                event = event[:-2]
                event += ";"
            else:
                event += operon + " " + str(absoluteIndex) + ";"
                genes.append(operon)
                indexes.append(absoluteIndex)
                absoluteIndex += 1
#    print event
    stopIndex = absoluteIndex

#    print "Performing inversions... Num before: %s Num after: %s" % (str(lengthBefore), str(lengthAfter))
#    print "Inversion:"
#    print formatGenome(before, after)
#    print ""

    eventType = "I"
    branchEvent = Event(eventType, indexes, genes, None, event)
    size = lengthBefore + lengthAfter + 1  #Additional +1 for terminus
    return branchEvent, indexes, genes, size

def reverseOperons(sequence):
    for index in sequence:
        if isinstance(index, list):
            index.reverse()

def performLoss(before, after, p):
    event = ""
    genes = []
    startPos = 0
    stopPos = 1

    if random.random() < 0.5:
        genome = before
        deleteBefore = True
    else:
        genome = after
        deleteBefore = False

    index = random.randint(0, len(genome)-1)
    if isinstance(genome[index], list):
        # print "is list"
        stopPos = min(geometricSampling(p), len(genome[index]))
        lossSection = genome[index][startPos:stopPos]

        if (stopPos - startPos) == len(genome[index]):
            del genome[index]
        else:
            del genome[index][startPos:stopPos]
            if len(genome[index]) == 0:
                del genome[index]
            elif len(genome[index]) == 1:
                newSingleton = genome[index][0]
                del genome[index]
                genome.insert(index, newSingleton)
    else:
        # print "is singleton"
        lossSection = genome[index:index+1]
        del genome[index]

    absoluteIndex = getAbsoluteIndex(deleteBefore, before, after, index, startPos)
    startIndex = absoluteIndex

    for gene in lossSection:
        event += gene + " " + str(absoluteIndex) + ","
        genes.append(gene)
        absoluteIndex += 1
#    print event
    stopIndex = absoluteIndex

    # event = "Performing deletion... From before: %s Index: %s  Section length: %s" % (str(deleteBefore), str(index), str(stopPos-startPos))
#    print "Performing deletion... From before: %s Index: %s  Section length: %s" % (str(deleteBefore), str(index), str(stopPos-startPos))
#    print formatGenome(before, after)
#    print ""

    eventType = "L"
    branchEvent = Event(eventType, range(startIndex, stopIndex), genes)
    return branchEvent, range(startIndex, stopIndex), genes

def performDuplication(before, after, p):
    event = ""
    genes = []
    startPos = 0
    stopPos = 0
    operonDup = False
    sectionReversed = False

    if random.random() < 0.5:
        genome = before
        dupFromBefore = True
    else:
        genome = after
        dupFromBefore = False

    index = random.randint(0, len(genome)-1)
    if isinstance(genome[index], list):
        # print "is list"
        operon = genome[index]
        stopPos = min(geometricSampling(p), len(operon))
        dupSection = operon[startPos:stopPos]

        if len(dupSection) == len(operon):
            operonDup = True
    else:
        # print "is singleton"
        dupSection = genome[index:index+1]

    beforeStartIndex = getAbsoluteIndex(dupFromBefore, before, after, index, startPos)
    beforeStopIndex = beforeStartIndex + len(dupSection)

    if random.random() < 0.5:
        genome = before
        targetBefore = True
    else:
        genome = after
        targetBefore = False

    if dupFromBefore and not targetBefore:
        dupSection.reverse()
        sectionReversed = True
        targetIndex, targetPos = getNormalPos(genome)
        # absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)
    elif not dupFromBefore and targetBefore:
        dupSection.reverse()
        sectionReversed = True
        targetIndex, targetPos = getNormalPos(genome)
        # absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)
    else:
        targetIndex, targetPos = getPosOutsideOfRange(genome, index, startPos, stopPos)
        # if targetBefore:
        #     # absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)

        # else:
        #     # absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)
    absoluteIndex = getAbsoluteIndex(targetBefore, before, after, targetIndex, targetPos)
    startIndex = absoluteIndex

    if isinstance(genome[targetIndex], list):
        # print "target in list"

        for gene in reversed(dupSection):
            # print "inserting " + gene
            genome[targetIndex].insert(targetPos, gene)

        # event = "".join([gene+" "+str(absoluteIndex)+"," for gene in dupSection])
        # print event
#        print "Performing duplication... From before: %s Index: %s  Section length: %s Target before: %s Target index: %s Target position: %s" % (str(dupFromBefore), str(index), str(len(dupSection)), str(targetBefore), str(targetIndex), str(targetPos))
    else:
        # print "target in singleton"
        if operonDup:
            genome.insert(targetIndex, dupSection)
        else:
            for gene in reversed(dupSection):
                # print "inserting " + gene
                genome.insert(targetIndex, gene)
#        print "Performing duplication... From before: %s Index: %s Section length: %s Target before: %s Target index: %s" % (str(dupFromBefore), str(index), str(len(dupSection)), str(targetBefore), str(targetIndex))

    for gene in dupSection:
        event += gene + " " + str(absoluteIndex) + ","
        genes.append(gene)
        absoluteIndex += 1
#    print event
    stopIndex = absoluteIndex

#    print "Duplication:"
#    print formatGenome(before, after)
#    print ""

    eventType = "D"
    branchEvent = Event(eventType, range(startIndex, stopIndex), genes, range(beforeStartIndex, beforeStopIndex))
    return branchEvent, range(startIndex, stopIndex), genes, range(beforeStartIndex, beforeStopIndex), sectionReversed

def getAbsoluteIndex(inBefore, before, after, targetIndex, targetPos):
    index  = 0

    if inBefore:
        for count1 in range(targetIndex):
            if isinstance(before[count1], list):
                index += len(before[count1])
                # print(before[count1])
                # print(len(before[count1]))
                # print "index after operon %d: %d" % (count1, index)
            else:
                index += 1
        index += targetPos + 1
    else:
        for item in before:
            if isinstance(item, list):
                index += len(item)
            else:
                index += 1
        # print "Length of before"
        # print str(index)

        for count2 in range(targetIndex):
            if isinstance(after[count2], list):
                index += len(after[count2])
                # print(after[count2])
                # print(len(after[count2]))
                # print "index after operon %d: %d" % (count2, index)
            else:
                index += 1
        index += targetPos + 2

    return index

def getPosOutsideOfRange(genome, index, start, stop):
    targetPos = 0

    while True:
        targetIndex = random.randint(0, len(genome)-1)

        if isinstance(genome[targetIndex], list):
            operon = genome[targetIndex]
            targetPos = random.randint(0, len(operon)-1)

            if targetIndex != index:
                break
            elif targetPos not in range(start, stop):
                break
        else:
            break

    return targetIndex, targetPos

def getNormalPos(genome):
    targetPos = 0
    targetIndex = random.randint(0, len(genome)-1)

    if isinstance(genome[targetIndex], list):
        operon = genome[targetIndex]
        targetPos = random.randint(0, len(operon)-1)

    return targetIndex, targetPos

def createAncestor(maxLength, numOperons):
    global beforeTerminus
    global afterTerminus
    beforeTerminus = []
    afterTerminus = []
    seqLength = 0
    currNumOperons = 0
    terminusAdded = False

    currentSequence = beforeTerminus
    while currNumOperons < numOperons or seqLength < maxLength:
        prob = random.random()

        if terminusAdded:
            rnaSequence = ["5S", "23S", "16S"]
        else:
            rnaSequence = ["16S", "23S", "5S"]

        if prob <= 0.6:
            if prob <= 0.07:
                operon = rnaSequence
            else:
                operon = createOperon()
                if prob <= 0.15:
                    if prob <= 0.11:
                        operon.extend(rnaSequence)
                    else:
                        operon = rnaSequence + operon
            if (seqLength + len(operon)) <= maxLength:
                seqLength += len(operon)
                currentSequence.append(operon)
                currNumOperons += 1
        else:
            if seqLength + 1 <= maxLength:
                currentSequence.append(random.choice(aminoAcids))
                seqLength += 1
        
        if seqLength == maxLength and currNumOperons < numOperons:
            beforeTerminus = []
            afterTerminus = []
            seqLength = 0
            currNumOperons = 0
            terminusAdded = False
            currentSequence = beforeTerminus

        if terminusAdded == False and seqLength >= maxLength/2:
            terminusAdded = True
            currentSequence = afterTerminus

def createOperon():
    operon = []
    size = max(geometricSampling(0.6), 2)

    for i in range(size):
        operon.append(random.choice(aminoAcids))

    return operon

def geometricSampling(p, minValue = 1):
    value = minValue
    while random.random() > p:
        value += 1
    return value

def createFile(fileName, newickTree):
#    print('Creating file %s...' % (fileName));    
    file = open(fileName, "w+")
    Phylo.write(newickTree,fileName, 'newick') #Write the tree to the output file
    file.close()
    
if __name__ == '__main__':
    #generateTests formatting
    #generateTests(OutputFolderName, treeFileName, numGenesForGenome, minNumOperons, numEventsPerBranch, dupProb, dupPValue, lossProb, loss_p, invProb, inv_p, subProb, transProb, trans_p)
    #OutputFolderName - folder where all your files will be output
    #treeFileName - tree you want to create
    #numGenesForGenome - number of genes each genome will have
    #minNumOperons - the genomes will have atleast this many operons
    #numEventsPerBranch - each branch will have exactly this many events
    
    #All probabilities determine how likely the event will occur. **Probabilities have to add up to 1
    #All pValues determine how many genes the event will affect. The higher the pValue, the lower the chance that the event will be large
    #Examples: code below creates genomes of size 25 with atleast 3 operons and 3 events for each genome. Only transpositions will occur though. None of the other events have a probability.
    generateTests("generatorTesting", "tree2Leaf.dnd", 25, 3, 3, 0.0, 0.7, 0.0, 0.7, 0.0, 0.7, 0.0, 1.0, 0.7)
     #Examples: code below creates genomes of size 100 with atleast 5 operons and 8 events for each genome. Only duplications and losses will occur.
    generateTests("generatorTesting", "tree2Leaf.dnd", 100, 5, 8, 0.5, 0.7, 0.5, 0.7, 0.0, 0.7, 0.0, 0.0, 0.7)
     #Examples: code below creates genomes of size 80 with atleast 4 operons and 6 events for each genome. All events have a chance of occuring.
    generateTests("generatorTesting", "tree2Leaf.dnd", 80, 4, 6, 0.4, 0.7, 0.4, 0.7, 0.05, 0.7, 0.05, 0.10, 0.7)