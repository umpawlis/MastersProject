from Bio import Phylo
import os.path
import multiset
import numpy as np
import xlsxwriter
import matplotlib.pyplot as plt

newickFileName = 'Anc27_subtree.dnd'
strains = []
ancestralCounter = 0
deletionCost = 1
substitutionCost = 2
codonCost = 0.5
trackingId = 0
duplicateOperonCounter = {}
deletionEventCounter = {}
duplicationEventCounter = {}

######################################################
# Tracking Event
# Parameters:
# Description: Stores information about a pair of orthologous operons
######################################################
class TrackingEvent(object):
    trackingEventId = 0
    score = -1
    genome1Name = ""
    genome2Name = ""
    genome1Operon = ""
    genome2Operon = ""
    genome1OperonIndex = -1
    genome2OperonIndex = -1
    ancestralOperon = ""
    technique = ""
    numLosses = 0
    operonEvents = None

    def __init__(self, trackingEventId, score, genome1Name, genome2Name, genome1Operon, genome2Operon, genome1OperonIndex, genome2OperonIndex, ancestralOperon, technique):
        self.trackingEventId = trackingEventId
        self.score = score
        self.genome1Name = genome1Name
        self.genome2Name = genome2Name
        self.genome1Operon = genome1Operon
        self.genome2Operon = genome2Operon
        self.genome1OperonIndex = genome1OperonIndex
        self.genome2OperonIndex = genome2OperonIndex
        self.ancestralOperon = ancestralOperon
        self.technique = technique
        self.lostEventIds = []

    def printTrackingEvent(self):
        if self.operonEvents != None:
            operonEventsString = self.operonEvents.toStringAlignmentResults()
            print("{ Tracking Event Id: %s, \nAlignment Score: %s, \nGenome 1: %s, \nGenome 2: %s, \nGenome 1 Operon: %s, \nGenome 2 Operon: %s, \nGenome 1 Operon Index: %s, \nGenome 2 Operon Index: %s, \nAncestral Operon: %s, \nTechnique: %s, \nLost Event Ids: %s, \nOperon Events: %s}"%(self.trackingEventId, self.score, self.genome1Name, self.genome2Name, self.genome1Operon, self.genome2Operon, self.genome1OperonIndex, self.genome2OperonIndex, self.ancestralOperon, self.technique, self.lostEventIds, operonEventsString))
        else:
            print("{ Tracking Event Id: %s, \nAlignment Score: %s, \nGenome 1: %s, \nGenome 2: %s, \nGenome 1 Operon: %s, \nGenome 2 Operon: %s, \nGenome 1 Operon Index: %s, \nGenome 2 Operon Index: %s, \nAncestral Operon: %s, \nTechnique: %s, \nLost Event Ids: %s}"%(self.trackingEventId, self.score, self.genome1Name, self.genome2Name, self.genome1Operon, self.genome2Operon, self.genome1OperonIndex, self.genome2OperonIndex, self.ancestralOperon, self.technique, self.lostEventIds))
    #####################################
    #Getters
    #####################################
    def getTrackingEventId(self):
        return self.trackingEventId
    def getScore(self):
        return self.score
    def getGenome1Name(self):
        return self.genome1Name
    def getGenome2Name(self):
        return self.genome2Name
    def getGenome1Operon(self):
        return self.genome1Operon
    def getGenome2Operon(self):
        return self.genome2Operon
    def getGenome1OperonIndex(self):
        return self.genome1OperonIndex
    def getGenome2OperonIndex(self):
        return self.genome2OperonIndex
    def getAncestralOperon(self):
        return self.ancestralOperon
    def getTechnique(self):
        return self.technique
    def getLostEventIds(self):
        return self.lostEventIds
    def getOperonEvents(self):
        return self.operonEvents

    #####################################
    #Setters
    #####################################
    def setTrackingEventId(self, trackingEventId):
        self.trackingEventId = trackingEventId
    def setScore(self, score):
        self.score = score
    def setGenome1Name(self, genome1Name):
        self.genome1Name = genome1Name
    def setGenome2Name(self, genome2Name):
        self.genome2Name = genome2Name
    def setGenome1Operon(self, genome1Operon):
        self.genome1Operon = genome1Operon
    def setGenome2Operon(self, genome2Operon):
        self.genome2Operon = genome2Operon
    def setGenome1OperonIndex(self, genome1OperonIndex):
        self.genome1OperonIndex = genome1OperonIndex
    def setGenome2OperonIndex(self, genome2OperonIndex):
        self.genome2OperonIndex = genome2OperonIndex
    def setAncestralOperon(self, ancestralOperon):
        self.ancestralOperon = ancestralOperon
    def setTechnique(self, technique):
        self.technique = technique
    def setLostEventIds(self, lostEventIds):
        self.lostEventIds = lostEventIds
    def setOperonEvents(self, operonEvents):
        self.operonEvents = operonEvents

######################################################
# Operon Events
# Parameters:
# Description: Stores information about the events that have occured in the operon
######################################################
class OperonEvents(object):
    numMatches = 0
    numCodonMismatches = 0
    numMismatches = 0
    numSubstitutions = 0
    operon1 = None
    operon2 = None
    matrix = None
    operon1GeneLosses = 0
    operon2GeneLosses = 0
    operon1GeneDuplicates = 0
    operon2GeneDuplicates = 0
    operon1Gaps = []
    operon2Gaps = []
    duplicateSizesOp1 = []
    duplicateSizesOp2 = []
    alignedGenesInOperon1 = []
    alignedGenesInOperon2 = []
    lossesDueToSlidingWindowMethodOperon1 = 0
    lossesDueToSlidingWindowMethodOperon2 = 0
    duplicationsDueToSlidingWindowMethodOperon1 = 0
    duplicationsDueToSlidingWindowMethodOperon2 = 0

    def __init__(self, numMatches, numCodonMismatches, numMismatches, numSubstitutions, operon1, operon2, matrix):
        self.numMatches = numMatches
        self.numCodonMismatches = numCodonMismatches
        self.numMismatches = numMismatches
        self.numSubstitutions = numSubstitutions
        self.operon1 = operon1
        self.operon2 = operon2
        self.matrix = matrix

    def toStringOperonEvents(self):
        return "(Num Matches = %s,\nNum Codon Mismatches = %s,\nNum Mismatches = %s,\nNum Substitutions = %s,\nOperon 1 = %s,\nOperon 2 = %s,\nOperon 1 Gene Losses: %s,\nOperon 2 Gene Losses: %s,\nOperon 1 Gene Duplications: %s,\nOperon 2 Gene Duplications: %s)" % (self.numMatches, self.numCodonMismatches, self.numMismatches, self.numSubstitutions, self.operon1, self.operon2, self.operon1GeneLosses, self.operon2GeneLosses, self.operon1GeneDuplicates, self.operon2GeneDuplicates)

    def toStringAlignmentResults(self):
        return "Results of Alignment: (Num Matches = %s,\nNum Codon Mismatches = %s,\nNum Mismatches = %s,\nNum Substitutions = %s,\nOperon 1 = %s,\nOperon 2 = %s)" % (self.numMatches, self.numCodonMismatches, self.numMismatches, self.numSubstitutions, self.operon1, self.operon2)

    #####Getters#####
    def getNumMatches(self):
        return self.numMatches
    def getNumCodonMismatches(self):
        return self.numCodonMismatches
    def getNumMismatches(self):
        return self.numMismatches
    def getNumSubstitutions(self):
        return self.numSubstitutions
    def getOperon1(self):
        return self.operon1
    def getOperon2(self):
        return self.operon2
    def getMatrix(self):
        return self.matrix
    def getOperon1GeneLosses(self):
        return self.operon1GeneLosses
    def getOperon2GeneLosses(self):
        return self.operon2GeneLosses
    def getOperon1GeneDuplicates(self):
        return self.operon1GeneDuplicates
    def getOperon2GeneDuplicates(self):
        return self.operon2GeneDuplicates
    def getOperon1Gaps(self):
        return self.operon1Gaps
    def getOperon2Gaps(self):
        return self.operon2Gaps
    def getDuplicateSizesOp1(self):
        return self.duplicateSizesOp1
    def getDuplicateSizesOp2(self):
        return self.duplicateSizesOp2
    def getAlignedGenesInOperon1(self):
        return self.alignedGenesInOperon1
    def getAlignedGenesInOperon2(self):
        return self.alignedGenesInOperon2
    def getLossesDueToSlidingWindowMethodOperon1(self):
        return self.lossesDueToSlidingWindowMethodOperon1
    def getLossesDueToSlidingWindowMethodOperon2(self):
        return self.lossesDueToSlidingWindowMethodOperon2
    def getDuplicationsDueToSlidingWindowMethodOperon1(self):
        return self.duplicationsDueToSlidingWindowMethodOperon1
    def getDuplicationsDueToSlidingWindowMethodOperon2(self):
        return self.duplicationsDueToSlidingWindowMethodOperon2

    #####Setters#####
    def setNumMatches(self, numMatches):
        self.numMatches = numMatches
    def setNumCodonMismatches(self, numCodonMismatches):
        self.numCodonMismatches = numCodonMismatches
    def setNumMismatches(self, numMismatches):
        self.numMismatches = numMismatches
    def setNumSubstitutions(self, numSubstitutions):
        self.numSubstitutions = numSubstitutions
    def setOperon1(self, operon1):
        self.operon1 = operon1
    def setOperon2(self, operon2):
        self.operon2 = operon2
    def setMatrix(self, matrix):
        self.matrix = matrix
    def setOperon1GeneLosses(self, operon1GeneLosses):
        self.operon1GeneLosses = operon1GeneLosses
    def setOperon2GeneLosses(self, operon2GeneLosses):
        self.operon2GeneLosses = operon2GeneLosses
    def setOperon1GeneDuplicates(self, operon1GeneDuplicates):
        self.operon1GeneDuplicates = operon1GeneDuplicates
    def setOperon2GeneDuplicates(self, operon2GeneDuplicates):
        self.operon2GeneDuplicates = operon2GeneDuplicates
    def setOperon1Gaps(self, operon1Gaps):
        self.operon1Gaps = operon1Gaps
    def setOperon2Gaps(self, operon2Gaps):
        self.operon2Gaps = operon2Gaps
    def setDuplicateSizesOp1(self, duplicateSizesOp1):
        self.duplicateSizesOp1 = duplicateSizesOp1
    def setDuplicateSizesOp2(self, duplicateSizesOp2):
        self.duplicateSizesOp2 = duplicateSizesOp2
    def setAlignedGenesInOperon1(self, alignedGenesInOperon1):
        self.alignedGenesInOperon1 = alignedGenesInOperon1
    def setAlignedGenesInOperon2(self, alignedGenesInOperon2):
        self.alignedGenesInOperon2 = alignedGenesInOperon2
    def setLossesDueToSlidingWindowMethodOperon1(self, lossesDueToSlidingWindowMethodOperon1):
        self.lossesDueToSlidingWindowMethodOperon1 = lossesDueToSlidingWindowMethodOperon1
    def setLossesDueToSlidingWindowMethodOperon2(self, lossesDueToSlidingWindowMethodOperon2):
        self.lossesDueToSlidingWindowMethodOperon2 = lossesDueToSlidingWindowMethodOperon2
    def setDuplicationsDueToSlidingWindowMethodOperon1(self, duplicationsDueToSlidingWindowMethodOperon1):
        self.duplicationsDueToSlidingWindowMethodOperon1 = duplicationsDueToSlidingWindowMethodOperon1
    def setDuplicationsDueToSlidingWindowMethodOperon2(self, duplicationsDueToSlidingWindowMethodOperon2):
        self.duplicationsDueToSlidingWindowMethodOperon2 = duplicationsDueToSlidingWindowMethodOperon2

######################################################
# Strain
# Parameters:
# Description: Stores information about a Strain from the newick tree
######################################################
class Strain(object):
    name = ""
    sequence = []
    genes = []
    descendants = []
    operonPositions = []
    singletonDict = {}
    trackingEvents = []
    hasData = False

    #Class constructor
    def __init__(self, name, sequence, descendants, operonPositions, singletonDict):
        self.name = name
        self.sequence = sequence
        self.descendants = descendants
        self.operonPositions = operonPositions
        self.singletonDict = singletonDict

    #Prints the Strain content
    def printStrain(self):
        print("{ Name: %s, Sequence: %s, Descendants: %s }" %(self.name, self.sequence, self.descendants))

    def getName(self):
        return self.name
    def getSequence(self):
        return self.sequence
    def getDescendants(self):
        return self.descendants
    def getGenes(self):
        return self.genes
    def setGenes(self, genes):
        self.genes = genes
    def getOperonPositions(self):
        return self.operonPositions
    def getSingletonDict(self):
        return self.singletonDict
    def getTrackingEvents(self):
        return self.trackingEvents
    def setTrackingEvents(self, trackingEvents):
        self.trackingEvents = trackingEvents
    def setHasData(self, hasData):
        self.hasData = hasData
    def getHasData(self):
        return self.hasData

######################################################
# processFileSequence
# Parameters: sequence - sequence to get the operons from
# Description: extracts a list of operons from a sequence
######################################################
def processFileSequence(sequence):
    geneList = []
    operonList = []
    singletonList = {}
    index = 0
    geneIndex = 0
    operonPositions = []

    while index < len(sequence):
        #Operon
        if (sequence[index] == '[') or (sequence[index] == '-' and sequence[index + 1] == '['):
            startIndex = index
            operonPositions.append(geneIndex)

            while sequence[index] != ']':
                if sequence[index] == ',':
                    geneIndex += 1
                index += 1

            #increment the index to include the ]
            index += 1
            operonList.append(sequence[startIndex:index])
            if sequence[startIndex] == '[':
                geneList.extend([gene.strip() for gene in sequence[startIndex+1:index-1].split(',')])
            else:
                geneList.extend([gene.strip() for gene in sequence[startIndex+2:index-1].split(',')])

        #Singleton
        if index < len(sequence) and sequence[index] == ',':
            geneIndex += 1
            if sequence[index+2] != '[' and sequence[index+2] != '<' and sequence[index+3] != '[':
                index += 1
                startIndex = index

                while sequence[index] != ',' and index < len(sequence)-1:
                    index += 1

                #Add singleton gene
                operonList.append(sequence[startIndex:index])
                operonPositions.append(geneIndex)

                singletonList[str(geneIndex)] = sequence[startIndex:index]
                geneList.append(sequence[startIndex:index].strip().replace("-", ""))
                index -= 1
        index += 1

    return operonList, operonPositions, singletonList, geneList

######################################################
# traverseNewickTree
# Parameters: node - The node that we want to process
# Description: Traverses a provided newick tree and reads in sequences if any (uses post traversal)
######################################################
def traverseNewickTree(node):
    #Global variables
    global strains
    global ancestralCounter
    #Local variables
    leftChildStrain = None
    rightChildStrain = None

    #Check if the clade has children
    if len(node.clades) > 0:
        leftChildStrain = traverseNewickTree(node.clades[0])
        if len(node.clades) > 1:
            rightChildStrain = traverseNewickTree(node.clades[1])

    #Check if the clade has a name, if it does, check if it has a directory for its sequence
    if node.name is not None and len(node.name) > 0:
        if os.path.isdir(node.name):
            if os.path.isfile(node.name + '/sequence.txt'):
                fileGeneSequence = open(node.name + '/sequence.txt', 'r').read()
                operons, operonPositions, singletonDict, allGenes = processFileSequence(fileGeneSequence)

                #Random newline characters in file are causing the alignment scores to be messed up
                for x in range(0, len(operons)):
                    operons[x] = (operons[x]).replace('\r', '').replace('\n', '')
                for x in range(0, len(allGenes)):
                    allGenes[x] = (allGenes[x]).replace('\r', '').replace('\n', '')
                for key in singletonDict.keys():
                    singletonDict[key] = (singletonDict[key]).replace('\r', '').replace('\n', '')

                strain = Strain(node.name, operons, [], operonPositions, singletonDict)
                strain.setGenes(allGenes)
                strain.setHasData(True)
                strains.append(strain)

                return strain

            else:
                print('No sequence file found for node: %s' % node.name)
        else:
            print('No directory found for node: %s' % node.name)

    if leftChildStrain is not None and leftChildStrain.getHasData() and rightChildStrain is not None and rightChildStrain.getHasData():
        print('These are the strains being compared: %s, %s'%(leftChildStrain.getName(), rightChildStrain.getName()))
        ancestralOperons, trackingEvents = processStrains(leftChildStrain, rightChildStrain)

        ancestralCounter += 1
        node.name = 'Ancestor %s' % (ancestralCounter)
        ancestor = Strain('Ancestor %s' % (ancestralCounter), ancestralOperons, [leftChildStrain.getName(), rightChildStrain.getName()], [], {})
        ancestor.setTrackingEvents(trackingEvents)
        ancestor.setHasData(False)
        strains.append(ancestor)

        return ancestor

    #If the left child has a sequence, return it
    elif leftChildStrain is not None and leftChildStrain.getHasData():
        return leftChildStrain

    #If the right child has a sequence, return it
    elif rightChildStrain is not None and rightChildStrain.getHasData():
        return rightChildStrain

    #If neither has a sequence, return None
    else:
        return None

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
# processStrains
# Parameters:
# Description:
######################################################
def processStrains(strain1, strain2):
    coverageTracker1 = {}
    coverageTracker2 = {}
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()

    for y in range(0, len(sequence1)):
        coverageTracker1[y] = False

    for x in range(0, len(sequence2)):
        coverageTracker2[x] = False

    #Global Alignment
    trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter = detectOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2)

    #Local Alignment
    localAlignmentTrackingEvents, coverageTracker1, coverageTracker2, localAlignmentCounter = detectOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2)
    if len(localAlignmentTrackingEvents) > 0:
        trackingEvents.append(localAlignmentTrackingEvents)

    #Handle singleton genes
    singletonDuplicatedG1, singletonLostG1, singletonTrackingEventsG1, coverageTracker1 = detectOrthologousSingletonGenes(coverageTracker1, strain1, strain1.getTrackingEvents())
    singletonDuplicatedG2, singletonLostG2, singletonTrackingEventsG2, coverageTracker2 = detectOrthologousSingletonGenes(coverageTracker2, strain2, strain2.getTrackingEvents())
    if len(singletonTrackingEventsG1) > 0:
        trackingEvents.append(singletonTrackingEventsG1)
    if len(singletonTrackingEventsG2) > 0:
        trackingEvents.append(singletonTrackingEventsG2)

    #Handle remaining operons
    operonDuplicateG1, operonLossG1, operonTrackingEventsG1, coverageTracker1 = detectDuplicateOperons(coverageTracker1, strain1)
    operonDuplicateG2, operonLossG2, operonTrackingEventsG2, coverageTracker2 = detectDuplicateOperons(coverageTracker2, strain2)
    if len(operonTrackingEventsG1) > 0:
        trackingEvents.append(operonTrackingEventsG1)
    if len(operonTrackingEventsG2) > 0:
        trackingEvents.append(operonTrackingEventsG2)

    #create dot plot
    if len(trackingEvents) > 0:
        createDotPlot(trackingEvents, strain1, strain2)
        updateGlobalTrackers(trackingEvents, strain1.getSequence(), strain2.getSequence())

    #Resolve remaining operons via self global alignment
    trackerDebugger(coverageTracker1, coverageTracker2, sequence1, sequence2)

    print('#' * 70)
    print('Statistics for the following strains: %s, %s' %(strain1.getName(), strain2.getName()))
    print('Total number of operons and singletons for %s: %s' %(strain1.getName(), len(coverageTracker1)))
    print('Total number of operons and singletons for %s: %s' %(strain2.getName(), len(coverageTracker2)))
    print('Number of orthologs found through global alignment: %s' %(globalAlignmentCounter))
    print('Number of orthologs found through local alignment: %s' %(localAlignmentCounter))
    print('Number of singletons identified as duplicates in %s: %s' %(strain1.getName(), singletonDuplicatedG1))
    print('Number of singletons identified as duplicates in %s: %s' %(strain2.getName(), singletonDuplicatedG2))
    print('Number of singletons lost in %s: %s' %(strain1.getName(), singletonLostG1))
    print('Number of singletons lost in %s: %s' %(strain2.getName(), singletonLostG1))
    print('Number of operons lost in %s: %s' %(strain1.getName(), operonLossG1))
    print('Number of operons lost in %s: %s' %(strain2.getName(), operonLossG1))
    print('Number of operons identified as duplicates in %s: %s' %(strain1.getName(), operonDuplicateG1))
    print('Number of operons identified as duplicates in %s: %s' %(strain2.getName(), operonDuplicateG2))
    print('#' * 70)

    return None, None

######################################################
# updateGlobalTrackers
# Parameters:
# Description:
######################################################
def updateGlobalTrackers(trackingEvents, sequence1, sequence2):

    formattedSequence1, operon1SequenceConversion = formatAllOperons(sequence1)
    formattedSequence2, operon2SequenceConversion = formatAllOperons(sequence2)

    for i in range(0, len(trackingEvents)):
        if trackingEvents[i].getTechnique() == '2 Genome Global Alignment':
            opEvents = trackingEvents[i].getOperonEvents()
            #Gets duplicate counts based on similar genes found within alignment
            if opEvents and opEvents.getDuplicateSizesOp1() and len(opEvents.getDuplicateSizesOp1()) > 0:
                for size in opEvents.getDuplicateSizesOp1():
                    updateDuplicationCounter(size)
            if opEvents and opEvents.getDuplicateSizesOp2() and len(opEvents.getDuplicateSizesOp2()) > 0:
                for size in opEvents.getDuplicateSizesOp2():
                    updateDuplicationCounter(size)

            #Computes losses and duplications for each set of gaps
            if opEvents and opEvents.getOperon1Gaps() and len(opEvents.getOperon1Gaps()) > 0:
                gaps = opEvents.getOperon1Gaps()
                totalLosses = 0
                for gap in gaps:
                    if len(gap) > 0:
                        #Compute duplications and losses along with their sizes
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(gap, formattedSequence1, trackingEvents[i].getGenome1OperonIndex())
                        totalLosses += numUniqueFound
                        #Update duplication counter
                        if len(duplicationSizes) > 0:
                            for size in duplicationSizes:
                                updateDuplicationCounter(size)
                        #Update loss counter
                        if len(deletionSizes) > 0:
                            for size in deletionSizes:
                                updateDeletionCounter(size)
                #end for
                opEvents.setLossesDueToSlidingWindowMethodOperon1(totalLosses)

            if opEvents and opEvents.getOperon2Gaps() and len(opEvents.getOperon2Gaps()) > 0:
                gaps = opEvents.getOperon2Gaps()
                totalLosses = 0
                for gap in gaps:
                    if len(gap) > 0:
                        #Uniques tells us number of losses
                        numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(gap, formattedSequence2, trackingEvents[i].getGenome2OperonIndex())
                        totalLosses += numUniqueFound
                        #Update duplication counter
                        if len(duplicationSizes) > 0:
                            for size in duplicationSizes:
                                updateDuplicationCounter(size)
                        #Update loss counter
                        if len(deletionSizes) > 0:
                            for size in deletionSizes:
                                updateDeletionCounter(size)
                #end for
                opEvents.setLossesDueToSlidingWindowMethodOperon2(totalLosses)

######################################################
# createDotPlot
# Parameters:
# Description:
######################################################
def createDotPlot(trackingEvents, strain1, strain2):
    #The green ones represent operons with no differences
    green_x_coord = []
    green_y_coord = []

    #Yellow ones represent scores between 1 and 2
    yellow_x_coord = []
    yellow_y_coord = []

    #Red ones represent scores between 3 and above
    red_x_coord = []
    red_y_coord = []

    #Blue ones represent a local alignment
    blue_x_coord = []
    blue_y_coord = []

    print("x" * 70)
    for i in range(0, len(trackingEvents)):
        if trackingEvents[i].getTechnique() == '2 Genome Global Alignment' or trackingEvents[i].getTechnique() == 'Local Alignment':
            #Assign the coords to the appropriate array
            if trackingEvents[i].getTechnique() == 'Local Alignment':
                blue_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                blue_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            elif trackingEvents[i].getScore() == 0:
                green_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                green_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            elif trackingEvents[i].getScore() == 1 or trackingEvents[i].getScore() == 2:
                yellow_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                yellow_y_coord.append(trackingEvents[i].getGenome2OperonIndex())
            else:
                red_x_coord.append(trackingEvents[i].getGenome1OperonIndex())
                red_y_coord.append(trackingEvents[i].getGenome2OperonIndex())

            print('x-axis: %s, y-axis: %s' %(trackingEvents[i].getGenome1OperonIndex(), trackingEvents[i].getGenome2OperonIndex()))

    #If we have any coordinates to plot, display them
    if len(green_x_coord) > 0 or len(yellow_x_coord) > 0 or len(red_x_coord) > 0 or len(blue_x_coord) > 0:
        f = plt.figure()
        plt.title("Orthologous Operon Mapping")
        plt.plot(green_x_coord, green_y_coord, 'go', yellow_x_coord, yellow_y_coord, 'yo', red_x_coord, red_y_coord, 'ro', blue_x_coord, blue_y_coord, 'bo')
        plt.axis([0, len(trackingEvents)+5, 0, len(trackingEvents)+5])
        plt.ylabel('Operon Position in Genome 1')
        plt.xlabel('Operon Position in Genome 2')
        plt.show()
        f.savefig("%s %s.pdf" %(strain1.getName(), strain2.getName()), bbox_inches='tight')
    else:
        print('No plot to display!')
    print("x" * 70)

######################################################
# detectDuplicateOperons
# Parameters:
# Description:
######################################################
def detectDuplicateOperons(coverageTracker, strain):
    global trackingId

    sequence = strain.getSequence()
    genomeName = strain.getName()
    operonDuplicate = 0
    operonLoss = 0
    numInvertedDuplicates = 0
    trackingEvents = []

    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False and len(sequence[i].split(',')) > 1:
            print('\n&&&&&&&&&& Duplicate Alignment &&&&&&&&&&&&&&&&')
            print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')

            duplicateEvent = duplicateAlignment(i, sequence[i], sequence, genomeName)
            coverageTracker[i] = True

            #checks if duplicate event is null
            if duplicateEvent:
                #Not null, we found a partner for this operon
                print('Found a matching duplicate operon, therefore not adding to ancestor!')
                duplicateEvent.printTrackingEvent()

                #Increment counter to indicate we successfully found a match
                operonDuplicate += 1

                #Check if inverted duplicates
                if ('-' in duplicateEvent.getGenome1Operon() and '-' not in duplicateEvent.getGenome2Operon()) or ('-' not in duplicateEvent.getGenome1Operon() and '-' in duplicateEvent.getGenome2Operon()):
                    numInvertedDuplicates += 1

            else:
                print('No duplicate ortholog found for operon: %s, therefore it will be added to ancestor as it is a loss' % (sequence[i]))
                trackingId += 1
                trackingEvent = TrackingEvent(trackingId, 0, genomeName, '', sequence[i], '', i, -1, sequence[i], "Duplicate Alignment (No match found)")
                trackingEvent = trackLossEvents(trackingEvent, strain.getTrackingEvents())

                #Indicates operon is a loss
                operonLoss += 1

                #decides whether to add the event or not
                if len(trackingEvent.getLostEventIds()) >= 2:
                    print('Removing this operon because it was lost two times in a row')
                else:
                    trackingEvents.append(trackingEvent)
                trackingEvent.printTrackingEvent()
            print('\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
            print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n')

    return operonDuplicate, operonLoss, trackingEvents, coverageTracker

######################################################
# duplicateAlignment
# Parameters:
# Description: Performs a global alignment on genome itself to find best source of duplicate
######################################################
def duplicateAlignment(g1OperonIndex, g1Operon, g1Sequence, genomeName1):
    print('Duplicate Alignment')

    optimalScore = -1
    duplicateEvent = None

    for x in range(0, len(g1Sequence)):
        #Ignore global alignment on operon itself
        if x != g1OperonIndex:
            #Compute the set differences
            setDifference, operon1, operon2, numDifferentGenes = formatAndComputeOperonDifferences(g1Operon, g1Sequence[x])

            #Checks if either operons at in the - orientation
            reverseOp1 = reverseSequence(g1Operon)
            reverseOp2 = reverseSequence(g1Sequence[x])

            #Reverse operons if needed to
            if reverseOp1 and reverseOp2 == False:
                operon1.reverse()
            if reverseOp2 and reverseOp1 == False:
                operon2.reverse()

            #Threshold to check if the sequences are worth comparing
            if setDifference <= (max(len(operon1), len(operon2))//3):
                lowestScore, operonEvents = performGlobalAlignment(operon1, operon2)

                if optimalScore == -1 or (lowestScore < optimalScore):
                    optimalScore = lowestScore
                    duplicateEvent = TrackingEvent(0, optimalScore, genomeName1, genomeName1, g1Operon, g1Sequence[x], g1OperonIndex, x, "", "Duplicate Alignment")
                    duplicateEvent.setOperonEvents(operonEvents)
                    distance = abs(g1OperonIndex - x)
                elif lowestScore == optimalScore and (abs(g1OperonIndex - x) < distance):
                    optimalScore = lowestScore
                    duplicateEvent = TrackingEvent(0, optimalScore, genomeName1, genomeName1, g1Operon, g1Sequence[x], g1OperonIndex, x, "", "Duplicate Alignment")
                    duplicateEvent.setOperonEvents(operonEvents)
                    distance = abs(g1OperonIndex - x)

    #check if we found a duplicate event
    if duplicateEvent:
        incrementDuplicateTracker(operon1)

    return duplicateEvent

######################################################
# incrementDuplicateTracker
# Parameters:
# Description: Increments the tracker according by using the operon length as the key
######################################################
def incrementDuplicateTracker(operon):
    key = str(len(operon)) #Figures out which key we need to increment

    if key in duplicateOperonCounter:
        duplicateOperonCounter[str(key)] = duplicateOperonCounter[str(key)] + 1
    else:
        duplicateOperonCounter[str(key)] = 1

######################################################
# detectOrthologousSingletonGenes
# Parameters:
# Description:
######################################################
def detectOrthologousSingletonGenes(coverageTracker, strain, descendantsTrackingEvents):
    global trackingId
    sequence = strain.getSequence()
    genomeName = strain.getName()
    singletonLost = 0
    singletonDuplicated = 0
    trackingEvents = []

    for i in range(0, len(coverageTracker)):
        if coverageTracker[i] == False and len(sequence[i].split(',')) == 1:
            addToAncestor, matchIndex, coverageTracker = resolveSingleton(sequence, i, coverageTracker)

            if addToAncestor:
                #If no match found then it's a loss so add to ancestor
                trackingId += 1
                trackingEvent = TrackingEvent(trackingId, 0, genomeName, '', sequence[i], '', i, -1, sequence[i], "Singleton Alignment")
                trackingEvent = trackLossEvents(trackingEvent, descendantsTrackingEvents)

                #Indicates a loss
                singletonLost += 1
                #decides whether to add the event or not
                if len(trackingEvent.getLostEventIds()) >= 2:
                    print('Removing this singleton because it was lost two times in a row')
                else:
                    trackingEvents.append(trackingEvent)
                trackingEvent.printTrackingEvent()
            else:
                #Increment Counter to indicate success in finding source
                singletonDuplicated += 1
    return singletonDuplicated, singletonLost, trackingEvents, coverageTracker

######################################################
# trackLossEvents
# Parameters:
# Description: Appends the previous tracking loss Ids to the current nodes tracking Id
######################################################
def trackLossEvents(currTrackingEvent, previousTrackingEvents):
    if len(previousTrackingEvents) == 0:
        #If there are no tracking events, then this is a leaf we're dealing with
        lostEventIds = currTrackingEvent.getLostEventIds()
        lostEventIds.append(currTrackingEvent.getTrackingEventId())
        currTrackingEvent.setLostEventIds(lostEventIds)
    else:
        #If there are tracking events, then this is not a leaf. Two cases to consider, it was lost in the previous tracking events or not
        for i in range(0, len(previousTrackingEvents)):
            if previousTrackingEvents[i].getAncestralOperon().strip() == currTrackingEvent.getAncestralOperon().strip():
                #We found the tracking event associated with the operon from the previous ancestral node
                lostEventIds = previousTrackingEvents[i].getLostEventIds()
                lostEventIds.append(currTrackingEvent.getTrackingEventId())
                currTrackingEvent.setLostEventIds(lostEventIds)
                #If there was a loss in the previous ancestor for the same operon, remove the operon from both the previous and current ancestor
                #if len(lostEventIds) >= 2:
                    #previousTrackingEvents.pop(i)

    return currTrackingEvent

######################################################
# resolveSingleton
# Parameters:
# Description: Finds original copy of singleton
######################################################
def resolveSingleton(sequence, singletonIndex, coverageTracker):
    addToAncestor = False
    sourceIndex = -1
    #Check if an exact gene exists in an operon
    for o in range(0, len(sequence)):
        #Don't compare to itself
        if o != singletonIndex:
            #Get a list of operon genes
            setDifference, singletonGene, operonGenes, numDifferentGenes = formatAndComputeOperonDifferences(sequence[singletonIndex], sequence[o])
            #If a match found and no other match found
            if (singletonGene[0] in operonGenes) and sourceIndex == -1:
                sourceIndex = o
                distance = abs(o - singletonIndex)
            #If a match found with smaller distance
            elif (singletonGene[0] in operonGenes) and distance > abs(o - singletonIndex):
                sourceIndex = o
                distance = abs(o - singletonIndex)

    coverageTracker[singletonIndex] = True
    if sourceIndex != -1:
        incrementDuplicateTracker(singletonGene)
        print('\n##### Singleton Source Found!! (NOT ADDED TO ANCESTOR) #####')
        print('The singleton gene %s was found in operon %s, index: %s' % (sequence[singletonIndex].strip(), sequence[sourceIndex].strip(), sourceIndex))
        print('###################################\n')
    else:
        addToAncestor = True
        print('\n##### No Singleton Source Found!!(ADDED TO ANCESTOR) #####')
        print('Could not find source for the singleton gene: %s' % (sequence[singletonIndex].strip()))
        print('###################################\n')

    return addToAncestor, sourceIndex, coverageTracker

######################################################
# updateDeletionCounter
# Parameters: deletionSize
# Description: Updates the deletion counter when a deletion event occurs.
######################################################
def updateDeletionCounter(deletionSize):
    global deletionEventCounter

    if str(deletionSize) in deletionEventCounter:
        deletionEventCounter[str(deletionSize)] += 1
    else:
        deletionEventCounter[str(deletionSize)] = 1

######################################################
# updateDuplicationCounter
# Parameters: duplicationSize
# Description: Updates the duplication counter when a duplication event occurs.
######################################################
def updateDuplicationCounter(duplicationSize):
    global duplicationEventCounter

    if str(duplicationSize) in duplicationEventCounter:
        duplicationEventCounter[str(duplicationSize)] += 1
    else:
        duplicationEventCounter[str(duplicationSize)] = 1

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
# detectOrthologsByLocalAlignment
# Parameters: Two descendants of the ancestor
# Description:
######################################################
def detectOrthologsByLocalAlignment(coverageTracker1, coverageTracker2, strain1, strain2):
    global trackingId

    #Finding optimal orthologs using local alignment
    minValue = min(len(coverageTracker1), len(coverageTracker2))
    trackingEvents = []
    localAlignmentCounter = 0
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()
    genomeName1 = strain1.getName()
    genomeName2 = strain2.getName()
    formattedSequence1, operonIndexes1 = formatAllOperons(strain1.getSequence())
    formattedSequence2, operonIndexes2 = formatAllOperons(strain2.getSequence())
    genesStrain1 = strain1.getGenes()
    genesStrain2 = strain2.getGenes()
    operonPositionList1 = strain1.getOperonPositions()
    operonPositionList2 = strain2.getOperonPositions()
    singletonDict1 = strain1.getSingletonDict()
    singletonDict2 = strain2.getSingletonDict()

    #Scan the matrix x times to find optimal local alignments everytime an entire matrix is scanned
    for x in range(0, minValue):
        highestScore = -1
        distance = 50   #arbitrary large number
        chosenOperon1 = []
        chosenOperon2 = []
        chosenStart = (0,0)
        chosenEnd = (0,0)
        chosenAligned1 = []
        chosenAligned2 = []

        for i in range(0, len(coverageTracker1)):
            #Check if operon was not covered and not a singleton
            if coverageTracker1[i] == False and len(sequence1[i].split(',')) > 1:
                for j in range(0, len(coverageTracker2)):
                    if coverageTracker1[i] == False and coverageTracker2[j] == False and len(sequence2[j].split(',')) > 1:
                        op1 = sequence1[i]
                        op2 = sequence2[j]

                        score, formattedOperon1, formattedOperon2, endPosition, startPosition, aligned1, aligned2, operonEvents = localAlignment(op1, op2, i, j, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2)

                        if (score > highestScore) or (score == highestScore and (abs(i - j)) < distance):
                            highestScore = score
                            rowIndex = i
                            colIndex = j
                            distance = abs(i - j)

                            chosenOperon1 = formattedOperon1
                            chosenOperon2 = formattedOperon2
                            chosenStart = startPosition
                            chosenEnd = endPosition
                            chosenAligned1 = aligned1
                            chosenAligned2 = aligned2
                            chosenOpEvents = operonEvents

        #After scanning the whole matrix if we found a best score, then store it
        if highestScore > -1:
            print('\n******** Local Alignment ***********')
            print('**************************************\n')

            localAlignmentCounter+=1
            coverageTracker1[rowIndex] = True
            coverageTracker2[colIndex] = True

            ancestralOperon = determineAncestor(chosenOperon1, chosenOperon2, chosenStart, chosenEnd, chosenAligned1, chosenAligned2, formattedSequence1, formattedSequence2, operonIndexes1[rowIndex], operonIndexes2[colIndex], chosenOpEvents)

            trackingId += 1
            trackingEvent = TrackingEvent(trackingId, highestScore, genomeName1, genomeName2, sequence1[rowIndex], sequence2[colIndex], rowIndex, colIndex, ancestralOperon, "Local Alignment")
            trackingEvent.printTrackingEvent()
            trackingEvents.append(trackingEvent)

            print('\n**************************************')
            print('**************************************\n\n')

    return trackingEvents, coverageTracker1, coverageTracker2, localAlignmentCounter

######################################################
# determineAncestor
# Parameters: op1, op2, startPosition, endPosition, aligned1, aligned2
# Description: Determines which ancestor to pick from local alignment results.
######################################################
def determineAncestor(op1, op2, startPosition, endPosition, aligned1, aligned2, sequence1, sequence2, opIndex1, opIndex2, operonEvents):
    chooseOp1 = True
    ancestor = []

    if len(op1) > len(op2):
        unaligned = getUnaligned(op1, startPosition[0]-1, endPosition[0]-1)
        #numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned2, unaligned)
        numUnique = len(op1)-len(aligned2)
        numDuplicateFound = 0
        alignedRange = (startPosition[0]-1, endPosition[0]-1)

        for geneList in unaligned:
            if geneList:
                numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(geneList, sequence1, opIndex1, alignedRange)
                numUnique -= numUniqueFound
                numDuplicateFound += numUniqueFound

        chooseOp1 = False

        if duplicationSizes:
            for size in duplicationSizes:
                updateDuplicationCounter(size)

        operonEvents.setOperon1GeneDuplicates(numDuplicateFound)
        operonEvents.setOperon2GeneLosses(numUnique)

        if deletionSizes:
            for size in deletionSizes:
                updateDeletionCounter(size)

    else:
        unaligned = getUnaligned(op2, startPosition[1]-1, endPosition[1]-1)
        # numUnique, unaligned, numDuplicateFound = compareDuplicates(aligned1, unaligned)
        numUnique = len(op2)-len(aligned1)
        numDuplicateFound = 0
        alignedRange = (startPosition[1]-1, endPosition[1]-1)

        for geneList in unaligned:
            if geneList:
                numUniqueFound, deletionSizes, duplicationSizes = findUniqueGenes(geneList, sequence2, opIndex2, alignedRange)
                numUnique -= numUniqueFound
                numDuplicateFound += numUniqueFound

        chooseOp1 = False
        operonEvents.setOperon1GeneLosses(numUnique)

        if deletionSizes:
            for size in deletionSizes:
                updateDeletionCounter(size)

        if duplicationSizes:
            for size in duplicationSizes:
                updateDuplicationCounter(size)

        operonEvents.setOperon2GeneDuplicates(numDuplicateFound)

    ancestor = unaligned[0] + aligned1 + unaligned[1]

    print(operonEvents.toStringOperonEvents())

    return ancestor

######################################################
# findUniqueGenes
# Parameters: geneList, sequence
# Description: Tries to find the list of genes in another operon of the genome.
######################################################
def findUniqueGenes(geneList, sequence, opIndex, alignedRange=(0,0)):
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

    while (comparisonSize >= 2) and genesNotFound:
        #Slide the window across the gene list, one gene at a time
        while (startIndex < comparisonSize) and (startIndex+comparisonSize <= len(geneList)) and genesNotFound:

            currentIndex = startIndex
            #Slide the window across the gene list by comparison size
            while currentIndex < len(geneList):
                if currentIndex+comparisonSize <= len(geneList):
                    gene = geneList[currentIndex:currentIndex+comparisonSize]
                    newSet = True
                elif not missedGene:
                    gene = geneList[currentIndex:len(geneList)]
                    newSet = True

                if newSet:
                    if geneInSequence(gene, sequence, len(gene), opIndex, alignedRange):
                        newRange = (currentIndex, currentIndex+comparisonSize-1)
                        if not checkOverlap(newRange, geneRanges):
                            geneRanges.append(newRange)
                            numGeneMatches += len(gene)
                            duplicationSizes.append(len(gene))
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
    if genesNotFound:
        inSet = True
        deletionSize = 0

        for index in reversed(range(len(geneList))):
            if not checkOverlap((index,index), geneRanges):
                if not inSet:
                    inSet = True
                deletionSize += 1
            else:
                geneList.pop(index)
                inSet = False
                if deletionSize != 0:
                    deletionSizes.append(deletionSize)
                    deletionSize = 0
        #Special case if all genes in list are losses
        if deletionSize != 0:
            deletionSizes.append(deletionSize)
            deletionSize = 0

    return numGeneMatches, deletionSizes, duplicationSizes

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
# getUnaligned
# Parameters: operon, startPosition, endPosition
# Description: Finds genes in the operon that are not part of the alignment
######################################################
def getUnaligned(operon, startPosition, endPosition):
    unaligned = []
    beforeList = []
    afterList = []

    for i in range(len(operon)):
        if i < startPosition:
            beforeList.append(operon[i])
        elif i > endPosition:
            afterList.append(operon[i])

    unaligned.append(beforeList)
    unaligned.append(afterList)

    return unaligned

######################################################
# localAlignment
# Parameters:
# Description: Performs local alignment on two operons
######################################################
def localAlignment(op1, op2, op1Position, op2Position, genesStrain1, genesStrain2, operonPositionList1, operonPositionList2, singletonDict1, singletonDict2):
    global fullAlignmentCounter
    #print('\nPerforming local alignment..')
    #print(op1)
    #print(op2)

    #Algorithm parameters
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5
    mismatch = -1.0
    gap = -1.0

    #Check if we need to reverse any of the operons
    reverseOp1 = reverseSequence(op1)
    reverseOp2 = reverseSequence(op2)

    #Compute the set differences between the two operons
    setDifference, operon1, operon2, numberOfDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

    leftAdjustment1 = 0
    leftAdjustment2 = 0
    rightAdjustment1 = 0
    rightAdjustment2 = 0
    #Reverse operons if needed to
    if reverseOp1:
        operon1.reverse()
    if reverseOp2:
        operon2.reverse()

    #Initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1))
    maxScore = 0
    extensionScore = 0
    returningScore = -1
    maxPosition = (0, 0)

    #Perform the Local Alignment
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):

            #check if genes are identical
            if operon1[a-1].split('_')[0].strip() == operon2[b-1].split('_')[0].strip():

                #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                if operon1[a-1].strip() == operon2[b-1].strip():
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + matchWithCodon
                else:
                    scoreMatrix[a][b] = scoreMatrix[a-1][b-1] + matchWithoutCodon
            else:
                scoreMatrix[a][b] = max(scoreMatrix[a-1][b] + gap, scoreMatrix[a][b-1] + gap, scoreMatrix[a-1][b-1] + mismatch, 0)

            if scoreMatrix[a][b] > maxScore:
                maxScore = scoreMatrix[a][b]
                maxPosition = (a, b)

    #Trace back score matrix to get alignment
    aligned1, aligned2, numGaps, endPosition, operonEvents = traceback(operon1, operon2, scoreMatrix, maxPosition)

    if len(operon1) <= len(operon2):
        shortestLength = len(operon1)
    else:
        shortestLength = len(operon2)

    #Adjust gene position we are looking at for extension based on the orienation of the operon
    if reverseOp1:
        leftAdjustment1 = len(operon1)-1
        leftAdjustment1 -= endPosition[0]-1
        rightAdjustment1 += (len(operon1) - (maxPosition[0]))
    else:
        leftAdjustment1 += endPosition[0]-1
        rightAdjustment1 = len(operon1)-1
        rightAdjustment1 -= (len(operon1) - (maxPosition[0]))

    if reverseOp2:
        leftAdjustment2 = len(operon2)-1
        leftAdjustment2 -= endPosition[1]-1
        rightAdjustment2 += (len(operon2) - (maxPosition[1]))
    else:
        leftAdjustment2 += endPosition[1]-1
        rightAdjustment2 = len(operon2)-1
        rightAdjustment2 -= (len(operon2) - (maxPosition[1]))

    #Only extend operons when the current alignment has no gaps
    if numGaps == 0 and genesStrain1 and genesStrain2:
        if len(aligned1) == shortestLength:
            fullAlignmentCounter += 1
            #print("One operon is a SUBSET of the other. Trying extention..")
            extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            returningScore = maxScore + extensionScore
            extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
            returningScore += extensionScore
            printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left and right extension")
        elif len(aligned1) < shortestLength:
            if maxPosition[0] == len(scoreMatrix)-1 or maxPosition[1] == len(scoreMatrix[0])-1:
                #print("Trying right extension..")
                extensionScore = extendAlignment("right", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], rightAdjustment1, rightAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after right extension")
                    returningScore = maxScore + extensionScore
            elif endPosition[0] == 1 or endPosition[1] == 1:
                #print("Trying left extension..")
                extensionScore = extendAlignment("left", operon1, operon2, genesStrain1, genesStrain2, operonPositionList1[op1Position], operonPositionList2[op2Position], leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2)
                if len(aligned1) >= shortestLength:
                    printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after left extension")
                    returningScore = maxScore + extensionScore
    #No extension but still qualifies because it has no gaps in alignment
    elif numGaps == 0 and len(aligned1) == shortestLength:
        fullAlignmentCounter += 1
        printAlignments(op1Position, op2Position, operon1, operon2, aligned1, aligned2, "after no extension due to missing gene list(s)")
        returningScore = maxScore
    #print("Final alignment:")
    #print(aligned1)
    #print(aligned2)

    return returningScore, operon1, operon2, maxPosition, endPosition, aligned1, aligned2, operonEvents

######################################################
# printStrains, printAlignments, printStats
# Parameters:
# Description: Prints a summary of the local alignments.
######################################################
def printStrains(strain1, strain2):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nResults of Local Alignment: %s & %s\n" %(strain1, strain2))
    file.write("...\n")
    file.close()

def printAlignments(op1Position, op2Position, operon1, operon2, alignment1, alignment2, message):
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("Strain 1 Operon %s: %s\n" %(op1Position, operon1))
    file.write("Strain 2 Operon  %s: %s\n" %(op2Position, operon2))
    file.write("Alignment Result (%s):\n" %message)
    file.write("%s\n" %alignment1)
    file.write("%s\n\n" %alignment2)

    file.close()

def printStats():
    global fullAlignmentCounter
    global extensionCounter
    file  = open('localAlignmentResults.txt', 'a+')

    file.write("\nTotals:\n")
    file.write("Number of full alignments that occurred: %s\n" %fullAlignmentCounter)
    file.write("Number of extensions performed: %s\n\n" %extensionCounter)

    file.close()
    fullAlignmentCounter = 0
    extensionCounter = 0

######################################################
# extendAlignment
# Parameters: operon1, operon2
# Description: Tries to extend the alignment of two operons by using their gene lists.
######################################################
def extendAlignment(direction, operon1, operon2, genesStrain1, genesStrain2, opGenePosition1, opGenePosition2, leftAdjustment1, leftAdjustment2, reverseOp1, reverseOp2, aligned1, aligned2, singletonDict1, singletonDict2):
    global extensionCounter
    mismatch = False
    operonRange1 = range(opGenePosition1, opGenePosition1+len(operon1))
    operonRange2 = range(opGenePosition2, opGenePosition2+len(operon2))
    opGenePosition1 += leftAdjustment1
    opGenePosition2 += leftAdjustment2

    extensionScore = 0
    matchWithCodon = 1.0
    matchWithoutCodon = 0.5

    while not mismatch:
        if direction == "left":
            if reverseOp1:
                opGenePosition1 += 1
            else:
                opGenePosition1 -= 1

            if reverseOp2:
                opGenePosition2 += 1
            else:
                opGenePosition2 -= 1
        elif direction == "right":
            if reverseOp1:
                opGenePosition1 -= 1
            else:
                opGenePosition1 += 1

            if reverseOp2:
                opGenePosition2 -= 1
            else:
                opGenePosition2 += 1

        if (opGenePosition1 in range(0, len(genesStrain1))) and (opGenePosition2 in range(0, len(genesStrain2))):
            #Check if the gene is a singlton or part of an operon
            if (str(opGenePosition1) in singletonDict1 and str(opGenePosition2) in singletonDict2):
                if (('-' in singletonDict1[str(opGenePosition1)]) != reverseOp1) or (('-' in singletonDict2[str(opGenePosition2)]) != reverseOp2):
                    mismatch = True
            elif (str(opGenePosition1) in singletonDict1 and opGenePosition2 in operonRange2):
                if (('-' in singletonDict1[str(opGenePosition1)]) != reverseOp1):
                    mismatch = True
            elif (str(opGenePosition2) in singletonDict2 and opGenePosition1 in operonRange1):
                if (('-' in singletonDict2[str(opGenePosition2)]) != reverseOp2):
                    mismatch = True
            else:
                mismatch = True

            if not mismatch:
                gene1 = genesStrain1[opGenePosition1]
                gene2  = genesStrain2[opGenePosition2]

                if gene1.split('_')[0].strip() == gene2.split('_')[0].strip():
                    #Codons match. Here we are comparing the genes with codons because if codons match, then whole gene will match
                    if gene1.strip() == gene2.strip():
                        extensionScore += matchWithCodon
                    else:
                        extensionScore += matchWithoutCodon

                    aligned1.insert(0, gene1)
                    aligned2.insert(0, gene2)
                    extensionCounter += 1
                else:
                    mismatch = True
        else:
            mismatch = True
    return extensionScore

######################################################
# traceback
# Parameters: operon1, operon2, scoreMatrix, startPosition
# Description: treaverses score matrix to determine optimal alignment
######################################################
def traceback(operon1, operon2, scoreMatrix, startPosition):

    match = 1
    codonMismatch = 0
    mismatch = 0
    substitution = 0
    operon1Losses = 0
    operon2Losses = 0
    operon1Duplications = 0
    operon2Duplications = 0

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = startPosition
    move = nextMove(scoreMatrix, x, y)
    numGaps = 0

    while move != END:
        if move == DIAG:
            if operon1[x-1] == operon2[y-1]:
                match += 1
            elif operon1[x-1].split('_')[0].strip() == operon2[y-1].split('_')[0].strip():
                codonMismatch += 1
            else:
                mismatch += 1

            aligned_seq1.append(operon1[x - 1])
            aligned_seq2.append(operon2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(operon1[x - 1])
            aligned_seq2.append('-')
            x -= 1
            numGaps += 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(operon2[y - 1])
            y -= 1
            numGaps += 1

        move = nextMove(scoreMatrix, x, y)

    aligned_seq1.append(operon1[x - 1])
    aligned_seq2.append(operon2[y - 1])
    endPosition = (x,y)

    operonEvents = OperonEvents(match, codonMismatch, mismatch, substitution, operon1, operon2, scoreMatrix, operon1Losses, operon2Losses, operon1Duplications, operon2Duplications)

    return list(reversed(aligned_seq1)), list(reversed(aligned_seq2)), numGaps, endPosition, operonEvents

######################################################
# nextMove
# Parameters: scoreMatrix, x, y
# Description: determines which direction to traverse score matrix
######################################################
def nextMove(scoreMatrix, x, y):
    diag = scoreMatrix[x - 1][y - 1]
    up   = scoreMatrix[x - 1][y]
    left = scoreMatrix[x][y - 1]

    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')

######################################################
# detectOrthologsByGlobalAlignment
# Parameters:
# Description:
######################################################
def detectOrthologsByGlobalAlignment(strain1, strain2, coverageTracker1, coverageTracker2):
    #Compute the global alignment matrix
    globalAlignmentMatrix, operonEventMatrix = computeGlobalAlignmentMatrix(strain1, strain2)
    trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter = scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, operonEventMatrix, coverageTracker1, coverageTracker2, strain1, strain2)

    return trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter

######################################################
# scanGlobalAlignmentMatrixForOrthologs
# Parameters:
# Description: Scans matrix and identifies orthologous operons
######################################################
def scanGlobalAlignmentMatrixForOrthologs(globalAlignmentMatrix, operonEventMatrix, coverageTracker1, coverageTracker2, strain1, strain2):
    global trackingId
    global codonCost

    genomeName1 = strain1.getName()
    genomeName2 = strain2.getName()
    sequence1 = strain1.getSequence()
    sequence2 = strain2.getSequence()
    maxValue = findMaxValueInMatrix(globalAlignmentMatrix)
    currentScoreSelected = 0
    globalAlignmentCounter = 0
    trackingEvents = []

    #Keep iterating util we find all the optimal scores (Finding orthologs using global alignment)
    while currentScoreSelected <= maxValue:
        #Prioritize the selection of operons with the same sign
        for i in range(0, len(globalAlignmentMatrix)):
            for j in range(0, len(globalAlignmentMatrix[i])):
                #Check if this is a * score, if both operons have not been marked off and if both are the same sign
                if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False) and (('-' in sequence1[i] and '-' in sequence2[j]) or ('-' not in sequence1[i] and '-' not in sequence2[j])):
                    score = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                    #Check if the score matches the scores we're currently looking for
                    if score == currentScoreSelected:
                        #We found an ortholog in the global alignment matrix
                        print('\n##### Global Alignment #####')
                        trackingId += 1
                        globalAlignmentCounter+=1
                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        if score == 0:
                            #We found a perfect match, doesn't matter which operon we pick
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, sequence1[i], "2 Genome Global Alignment")
                        else:
                            #We found orthologs that are not a perfect match and need to be resolved
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, '', "2 Genome Global Alignment")
                        #Add the event to the tracking events list
                        trackingEvent.setOperonEvents(operonEventMatrix[i][j])
                        trackingEvents.append(trackingEvent)
                        trackingEvent.printTrackingEvent()
                        print('###################################\n')

        #Select the remaining operons with the optimal score
        for i in range(0, len(globalAlignmentMatrix)):
            for j in range(0, len(globalAlignmentMatrix[i])):
                #Check if this is a * score and if both operons have not been marked off
                if ('*' in str(globalAlignmentMatrix[i][j])) and (coverageTracker1[i] == False) and (coverageTracker2[j] == False):
                    score = float(str(globalAlignmentMatrix[i][j]).replace('*', ''))
                    #Check if the score matches the scores we're currently looking for
                    if score == currentScoreSelected:
                        #We found an ortholog in the global alignment matrix
                        print('\n##### Global Alignment #####')
                        trackingId += 1
                        globalAlignmentCounter+=1
                        coverageTracker1[i] = True
                        coverageTracker2[j] = True

                        if score == 0:
                            #We found a perfect match, doesn't matter which operon we pick
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, sequence1[i], "2 Genome Global Alignment")
                        else:
                            trackingEvent = TrackingEvent(trackingId, score, genomeName1, genomeName2, sequence1[i], sequence2[j], i, j, '', "2 Genome Global Alignment")
                        #Add the event to the tracking events list
                        trackingEvent.setOperonEvents(operonEventMatrix[i][j])
                        trackingEvents.append(trackingEvent)
                        trackingEvent.printTrackingEvent()
                        print('###################################\n')
        currentScoreSelected += codonCost

    return trackingEvents, coverageTracker1, coverageTracker2, globalAlignmentCounter

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

            #compute the set differences between the two operons
            setDifference, operon1, operon2, numDifferentGenes = formatAndComputeOperonDifferences(op1, op2)

            #Checks if either operons at in the - orientation
            reverseOp1 = reverseSequence(op1)
            reverseOp2 = reverseSequence(op2)

            #Reverse operons if needed to
            if reverseOp1 and reverseOp2 == False:
                operon1.reverse()
            if reverseOp2 and reverseOp1 == False:
                operon2.reverse()

            #Case 1: We have two singleton operons
            if len(operon1) == 1 and len(operon2) == 1:
                #Perfect match
                if operon1 == operon2:
                    globalAlignmentMatrix[x][y] =  str(0) + '*'
                #Mismatch
                else:
                    globalAlignmentMatrix[x][y] = -999

            #Case 2: Only one of them is a singleton operon
            elif (len(operon1) == 1 and len(operon2) > 1) or (len(operon2) == 1 and len(operon1) > 1):
                globalAlignmentMatrix[x][y] = -999

            #Case 3: None of them are singleton operons, perform a global alignment
            elif len(op1) > 1 and len(op2) > 1:
                score, operonEvents = performGlobalAlignment(operon1, operon2)

                globalAlignmentMatrix[x][y] = score
                operonEventMatrix[x][y] = operonEvents

                threshold = max(len(operon1), len(operon2))
                threshold = threshold//3

                if setDifference <= threshold:
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

######################################################
# outputResultsToExcel
# Parameters: strain1, strain2, sequence1, sequence2, resultMatrix
# Description: outputs the results into an excel file
######################################################
def outputResultsToExcel(strain1, strain2, sequence1, sequence2, resultMatrix):

    rowIndex = 1

    #Intialize the workbook
    workbook = xlsxwriter.Workbook('excel_files/Comparison - %s & %s.xlsx' % (strain1, strain2))
    worksheet = workbook.add_worksheet()
    worksheet2 = workbook.add_worksheet()

    #Used to emphasize potentially interesting operons
    cellFormat = workbook.add_format()
    cellFormat.set_pattern(1)
    cellFormat.set_bg_color('lime')

    #Write strain identifiers to the excel file
    worksheet.write(5, 0, strain1)
    worksheet.write(0, 5, strain2)
    worksheet2.write(5, 0, strain1)
    worksheet2.write(0, 5, strain2)

    #Write the operons to the excel file
    for x in range (0, len(sequence1)):
        worksheet.write((x + 2), 1, sequence1[x])
        worksheet2.write((x + 2), 1, sequence1[x])

    for x in range (0, len(sequence2)):
        worksheet.write(1, (x + 2), sequence2[x])
        worksheet2.write(1, (x + 2), sequence2[x])

    #Write the data from the matrix to the excel file
    for row in resultMatrix:

        #Track the excel indexes
        rowIndex += 1
        colIndex = 2

        for value in row:
            if '*' in str(value):
                worksheet.write(rowIndex, colIndex, str(value), cellFormat)
            else:
                worksheet.write(rowIndex, colIndex, str(value))
            colIndex += 1

    #Close the excel file
    workbook.close()

######################################################
# reverseSequence
# Parameters: operon - operon to check
# Description: checks if the operon needs to be reversed
######################################################
def reverseSequence(operon):

    reverseOperon = False

    if '-' in operon:
        reverseOperon = True

    return reverseOperon

######################################################
# formatAndComputeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: computes the set difference between two operons
# and processes the operon lists
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
# trackerDebugger
# Parameters:
# Description: Prints the content of each respective tracker
######################################################
def trackerDebugger(coverageTracker1, coverageTracker2, sequence1, sequence2):
    print('Remaining operons from each respective tracker:')
    for x in range(0, len(coverageTracker1)):
        if coverageTracker1[x] == False:
            print ('Sequence 1, index: %s, Operon: %s' % (x, sequence1[x]))
    for x in range (0, len(coverageTracker2)):
        if coverageTracker2[x] == False:
            print('Sequence 2, index: %s, Operon: %s' % (x, sequence2[x]))
    print('Finished printing trackers\n')

######################################################
# performGlobalAlignment
# Parameters:
# Description: Performs a forward and reversed global alignment and returns the best score
######################################################
def performGlobalAlignment(operon1, operon2):
    global substitutionCost
    global deletionCost
    global codonCost

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
                    scoreMatrix[a][b] = min(scoreMatrix[a-1][b-1] + codonCost, scoreMatrix[a-1][b] + deletionCost, scoreMatrix[a][b-1] + deletionCost, scoreMatrix[a-1][b-1] + substitutionCost)
            else:
                scoreMatrix[a][b] = min(scoreMatrix[a-1][b] + deletionCost, scoreMatrix[a][b-1] + deletionCost, scoreMatrix[a-1][b-1] + substitutionCost)

    #Compute the number of events that occured between the operons
    operonEvents = globalAlignmentTraceback(scoreMatrix, operon1, operon2)

    return scoreMatrix[len(operon1)][len(operon2)], operonEvents

######################################################
# globalAlignmentTraceback
# Parameters:
# Description: Performs a traceback on a given matrix
######################################################
def globalAlignmentTraceback(matrix, operon1, operon2):
    global substitutionCost
    global deletionCost
    global codonCost

    i = len(operon1)
    j = len(operon2)

    match = 0
    codonMismatch = 0
    mismatch = 0
    substitution = 0

    operon1Gaps = []
    operon1Gap = []
    operon1ConsecutiveGap = False #Tracks consecutive gaps

    operon2Gaps = []
    operon2Gap = []
    operon2ConsecutiveGap = False #Tracks consecutive gaps

    #Track the alignment
    alignmentSequence1 = []
    alignmentSequence2 = []

    while i > 0 or j > 0:
        #Case 1: Perfect match
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] and operon1[i-1] == operon2[j-1]:
            match += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + codonCost) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            codonMismatch += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 3: Substitution
        elif i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + substitutionCost):
            substitution += 1
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False
            alignmentSequence1.append(operon1[i-1])
            alignmentSequence2.append(operon2[j-1])
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and matrix[i][j] == (matrix[i-1][j] + deletionCost):
            index = i-1
            mismatch += 1
            i -= 1
            operon1ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon2ConsecutiveGap:
                operon2Gap.append(operon1[index])
                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.append(operon2Gap)
                operon2Gap = []
                operon2Gap.append(operon1[index])
                operon2ConsecutiveGap = True

        #Case 5: Mismatch - Gap in operon 1
        else:
            index = j - 1
            mismatch += 1
            j -= 1
            operon2ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon1ConsecutiveGap:
                operon1Gap.append(operon2[index])
                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.append(operon1Gap)
                operon1Gap = []
                operon1Gap.append(operon2[index])
                operon1ConsecutiveGap = True

    #Empty any remaining gaps
    if len(operon1Gap) > 0:
        operon1Gaps.append(operon1Gap)
        operon1Gap = []
    if len(operon2Gap) > 0:
        operon2Gaps.append(operon2Gap)
        operon2Gap = []

    #Computes number of unique genes and if the genes are not unique then they are removed from the gap
    operon1Gaps, geneDuplicateSizesInOperon1 = checkForMatchesInAlignment(operon1Gaps, alignmentSequence1)
    operon2Gaps, geneDuplicateSizesInOperon2 = checkForMatchesInAlignment(operon2Gaps, alignmentSequence2)

    operonEvents = OperonEvents(match, codonMismatch, mismatch, substitution, operon1, operon2, matrix)

    #Sets the information we need to perform the sliding window later
    operonEvents.setAlignedGenesInOperon1(alignmentSequence1)
    operonEvents.setAlignedGenesInOperon2(alignmentSequence2)
    operonEvents.setOperon1Gaps(operon1Gaps)
    operonEvents.setOperon2Gaps(operon2Gaps)
    operonEvents.setDuplicateSizesOp1(geneDuplicateSizesInOperon1)
    operonEvents.setDuplicateSizesOp2(geneDuplicateSizesInOperon2)

    return operonEvents

######################################################
# checkForMatchesInAlignment
# Parameters:
# Description: Takes an array of gaps and an alignment, then checks if the genes in the gap match any of the genes in the alignment, if they do then genes are popped off the gap array
######################################################
def checkForMatchesInAlignment(arrayOfGaps, alignedGenes):
    geneDuplicateSizes = []

    for gap in arrayOfGaps:
        #Initialize Window
        windowSize = len(gap)
        startIndex = 0
        endIndex = len(gap)

        #print('Current gap %s' % (gap))
        while windowSize > 1:
            genes = gap[startIndex:endIndex]

            #print('Current Window %s' %(genes))
            genesMatched = 0

            for x in range(0, len(alignedGenes)):
                if len(genes) > 0 and genes[0] == alignedGenes[x] and genesMatched == 0:
                    for y in range(0, len(genes)):
                        if (x+y) < len(alignedGenes) and genes[y] == alignedGenes[x+y]:
                            genesMatched +=1
                    if genesMatched != len(genes):
                        genesMatched = 0

            if genesMatched != 0 and genesMatched == len(genes):
                #print("Duplicate")
                #updateDuplicationCounter(len(genes))
                geneDuplicateSizes.append(len(genes))
                del gap[startIndex:endIndex]
                startIndex = endIndex
            else:
                startIndex+=1

            if (startIndex + windowSize) > len(gap):
                #reduce and reset
                windowSize = min(windowSize-1, len(gap))
                startIndex = 0
            endIndex = startIndex + windowSize

    return arrayOfGaps, geneDuplicateSizes

######################################################
# drawDuplicationLossDistributionPlot
# Parameters:
# Description: Constructs the Duplication and Loss Distribution Plot
######################################################
def drawDuplicationLossDistributionPlot(duplicateOperonCounter, duplicationEventCounter, deletionEventCounter):
    fig = plt.figure()
    w = 0.1
    print('Constructing Duplication and Loss Distribution Plot...')
    if len(duplicateOperonCounter) > 0:
        print('Total Duplicate Operon Results:')
        operonDuplication_x_coords, operonDuplication_y_coords = getCoordinates(duplicateOperonCounter)
        operonDuplication_x_coords = np.asarray(operonDuplication_x_coords)
        plt.bar(operonDuplication_x_coords - w, operonDuplication_y_coords, width=w, color='#e74c3c', align='center')

    if len(duplicationEventCounter) > 0:
        print('Total Duplicate Gene Results:')
        geneDuplication_x_coords, geneDuplication_y_coords = getCoordinates(duplicationEventCounter)
        geneDuplication_x_coords = np.asarray(geneDuplication_x_coords)
        plt.bar(geneDuplication_x_coords, geneDuplication_y_coords, width=w, color='#f1c40f', align='center')

    if len(deletionEventCounter) > 0:
        print('Total Results of Loss Gene Tracker:')
        geneLoss_x_coords, geneLoss_y_coords = getCoordinates(deletionEventCounter)
        geneLoss_x_coords = np.asarray(geneLoss_x_coords)
        plt.bar(geneLoss_x_coords + w, geneLoss_y_coords, width=w, color='#3498db', align='center')

    #Formats the axis
    plt.ylabel('Number of Events')
    plt.xlabel('Size of Sequence')
    plt.title('Destribution of Duplications and Losses')
    fig.set_size_inches(5, 8, forward=True)
    plt.show()
    fig.savefig("Duplicate_Tracker.pdf", bbox_inches='tight')
    print('Done Constructing Duplication and Loss Distribution Plot')


######################################################
# getCoordinates
# Parameters: mapper
# Description: iterates through a dictionary and gathers the x and y coordinates
######################################################
def getCoordinates(mapper):
    x_coords = []
    y_coords = []

    key = 1
    itemsIterated = 0
    while itemsIterated < len(mapper):
        if mapper.get(str(key)):
            value =  mapper.get(str(key))
            x_coords.append(key)
            y_coords.append(value)
            print("Size: %s => Num Events: %s" % (key, value))
            itemsIterated+=1
        else:
            x_coords.append(key)
            y_coords.append(0)
        key+=1
    return x_coords, y_coords

######################################################
#                       main
######################################################
print('Reading in newick tree from file: %s...' % (newickFileName))
newickTree = Phylo.read(newickFileName, 'newick')
Phylo.draw(newickTree)

#Traverses the newick tree to reconstruct the ancestral genomes
result = traverseNewickTree(newickTree.clade)

#Construct the distribution loss bar graph
global duplicateOperonCounter
global duplicationEventCounter
global deletionEventCounter
drawDuplicationLossDistributionPlot(duplicateOperonCounter, duplicationEventCounter, deletionEventCounter)
print('Ending of processing')