import globals
import copy
import multiset
import numpy as np
import matplotlib.pyplot as plt

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
# updateGlobalCodonMismatchCounter
# Parameters: strain
# Description: Increments the global counter for codon mismatches
######################################################
def updateGlobalCodonMismatchCounter(strain):
    tempString = copy.deepcopy(strain.codonMismatchDetails)
    tempString = tempString.replace('Codon Mismatch:', '').strip()
    if tempString:
        array = filter(None, tempString.split(';'))
        if len(array) > 0:
            globals.codonMismatchCounter += len(array)
        
######################################################
# updateGlobalCodonMismatchCounter
# Parameters: strain
# Description: Increments the global counter for codon mismatches
######################################################
def updateGlobalSubstitutionCounter(strain):
    tempString = copy.deepcopy(strain.substitutionDetails)
    tempString = tempString.replace('Substitution:', '').strip()
    if tempString:
        array = filter(None, tempString.split(';'))
        if len(array) > 0:
            globals.substitutionCounter += len(array)
            
######################################################
# updateGlobalInversionSizeDistributionCounter
# Parameters: strain
# Description: Increments the global size distributions for inversions
######################################################
def updateGlobalInversionSizeDistributionCounter(strain):
    if len(strain.inversionCounts) > 0:
        for size, count in strain.inversionCounts.items():
            if size in globals.inversionSizeDistributionCounter:
                globals.inversionSizeDistributionCounter[size] += count
            else:
                globals.inversionSizeDistributionCounter[size] = count
                
######################################################
# updateGlobalTranspositionSizeDistributionCounter
# Parameters: strain
# Description: Increments the global size distributions for transpositions
######################################################
def updateGlobalTranspositionSizeDistributionCounter(strain):
    if len(strain.transpositionCounts) > 0:
        for size, count in strain.transpositionCounts.items():
            if size in globals.transpositionSizeDistributionCounter:
                globals.transpositionSizeDistributionCounter[size] += count
            else:
                globals.transpositionSizeDistributionCounter[size] = count

######################################################
# updateGlobalInvertedTranspositionSizeDistributionCounter
# Parameters: strain
# Description: Increments the global size distributions for inverted transpositions
######################################################
def updateGlobalInvertedTranspositionSizeDistributionCounter(strain):
    if len(strain.invertedTranspositionCounts) > 0:
        for size, count in strain.invertedTranspositionCounts.items():
            if size in globals.invertedTranspositionSizeDistributionCounter:
                globals.invertedTranspositionSizeDistributionCounter[size] += count
            else:
                globals.invertedTranspositionSizeDistributionCounter[size] = count

######################################################
# updateGlobalDeletionCounter
# Parameters: strain
# Description: Increments the global deletion counter based on the strains deletion sizes
######################################################
def updateGlobalDeletionCounter(strain):
    if len(strain.deletionCounts) > 0:
        for size, count in strain.deletionCounts.items():
            if size in globals.deletionSizeCounter:
                globals.deletionSizeCounter[size] += count
            else:
                globals.deletionSizeCounter[size] = count

######################################################
# updateGlobalDuplicationCounter
# Parameters: strain
# Description: Increments the global duplication counter based on the strains duplication sizes
######################################################
def updateGlobalDuplicationCounter(strain):
    if len(strain.duplicationCounts) > 0:
        for size, count in strain.duplicationCounts.items():
            if size in globals.duplicationSizeCounter:
                globals.duplicationSizeCounter[size] += count
            else:
                globals.duplicationSizeCounter[size] = count

######################################################
# addDuplicationEventsToStrain
# Parameters: stain, duplication sizes, details about the description of the event
# Description: Increments the duplication counter and adds the details to the strain
######################################################
def addDuplicationEventsToStrain(strain, duplicationSizes, duplicationDescription):

    if duplicationDescription != None and len(duplicationSizes) > 0:
        strain.duplicationDetails += duplicationDescription
        for x in range(0, len(duplicationSizes)):
            if duplicationSizes[x] in strain.duplicationCounts:
                strain.duplicationCounts[duplicationSizes[x]] += 1
            else:
                strain.duplicationCounts[duplicationSizes[x]] = 1
    return strain

######################################################
# addDeletionEventsToStrain
# Parameters: stain, deletion sizes, details about the description of the event
# Description: Increments the deletion counter and adds the details to the strain
######################################################
def addDeletionEventsToStrain(strain, deletionSizes, deletionDescription):

    if deletionDescription != None and len(deletionSizes) > 0:
        strain.deletionDetails += deletionDescription
        for x in range(0, len(deletionSizes)):
            if deletionSizes[x] in strain.deletionCounts:
                strain.deletionCounts[deletionSizes[x]] += 1
            else:
                strain.deletionCounts[deletionSizes[x]] = 1
    return strain

######################################################
# createDotPlot
# Parameters:
# Description: Constructs a dot plot to indicate a mapping of the orthologous operons
######################################################
def createDotPlot(events, strain1, strain2, testFolder = ''):
    #Stores all of the coordinates
    x_coord = []
    y_coord = []

    #The black ones represent the origin and terminus
    black_x_coord = []
    black_y_coord = []

    #The green ones represent operons with no differences (Global alignment)
    green_x_coord = []
    green_y_coord = []

    #Yellow ones represent scores between 1 and 2 (Global alignment)
    yellow_x_coord = []
    yellow_y_coord = []

    #Orange ones represent scores between 3 and above (Global alignment)
    orange_x_coord = []
    orange_y_coord = []

    #Red ones represent a local alignment
    red_x_coord = []
    red_y_coord = []
    if globals.printToConsole:
        print("x" * 70)
    for i in range(0, len(events)):
        if events[i].technique == 'Global Alignment' or events[i].technique == 'Local Alignment':
            #Assign the coords to the appropriate array, the index represents the position of the operon with respect to the genome
            totalNumEvents = events[i].numMismatches + events[i].numCodonMismatches + events[i].numSubstitutions
            if events[i].technique == 'Local Alignment':
                red_x_coord.append(events[i].fragmentDetails1.point)
                red_y_coord.append(events[i].fragmentDetails2.point)
            elif events[i].fragmentDetails1.description == 'Terminus' or events[i].fragmentDetails1.description == 'Origin':
                black_x_coord.append(events[i].fragmentDetails1.point)
                black_y_coord.append(events[i].fragmentDetails2.point)
            elif totalNumEvents == 0:
                green_x_coord.append(events[i].fragmentDetails1.point)
                green_y_coord.append(events[i].fragmentDetails2.point)
            elif totalNumEvents == 1 or totalNumEvents == 2:
                yellow_x_coord.append(events[i].fragmentDetails1.point)
                yellow_y_coord.append(events[i].fragmentDetails2.point)
            else:
                orange_x_coord.append(events[i].fragmentDetails1.point)
                orange_y_coord.append(events[i].fragmentDetails2.point)

            #Get all coordinates into a single array
            x_coord.append(events[i].fragmentDetails1.point)
            y_coord.append(events[i].fragmentDetails2.point)
            #print('x-axis: %s, y-axis: %s' %(events[i].fragmentDetails1.point, events[i].fragmentDetails2.point))
    if len(black_x_coord) != 2:
        if globals.printToConsole:
            print('BREAK!')
    
    #If we have any coordinates to plot, display them
    if len(green_x_coord) > 0 or len(yellow_x_coord) > 0 or len(orange_x_coord) > 0 or len(red_x_coord) > 0:
        f = plt.figure()
        plt.title("Orthologous Operon Mapping")
        plt.plot(green_x_coord, green_y_coord, 'o', color = 'green')
        plt.plot( yellow_x_coord, yellow_y_coord, 'o', color = 'gold')
        plt.plot(orange_x_coord, orange_y_coord, 'o', color = 'orange')
        plt.plot(red_x_coord, red_y_coord, 'o', color = 'red')
        plt.plot(black_x_coord, black_y_coord, 'o', color = 'black')
        plt.axis([0, len(events)+5, 0, len(events)+5])
        plt.ylabel('Operon Position in %s' % (strain1.name))
        plt.xlabel('Operon Position in %s' % (strain2.name))
        plt.show()
        f.savefig(testFolder + '%s %s.pdf' %(strain1.name, strain2.name), bbox_inches='tight')
    else:
        if globals.printToConsole:
            print('No plot to display!')
    if globals.printToConsole:
        print("x" * 70)


######################################################
# normalizeIndexesForDotPlot
# Parameters:
# Description: Adjusts the dot plot points according to number of duplications and losses
######################################################
def normalizeIndexesForDotPlot(events, dup1, dup2, strain1, strain2):
    points = []
    lostPoints = []

    #Split the events into two lists, a points list which will be displayed in the dot plot and a lost points list for the lost operons
    for event in events:
        if event.score == -1:
            lostPoints.append(event)
        else:
            points.append(event)

    #Normalize the points in the points list
    points.sort(key=lambda x:x.fragmentDetails1.fragmentIndex, reverse=False)
    for x in range(0, len(points)):
        point = points[x]
        
        #These points are different, need to adjust it according to the losses and duplications
        index1 = point.fragmentDetails1.fragmentIndex
        index2 = point.fragmentDetails2.fragmentIndex

        lossCount1, dupCount1 = CountNumLossesAndDuplications(index1, lostPoints, dup1, strain1)
        lossCount2, dupCount2 = CountNumLossesAndDuplications(index2, lostPoints, dup2, strain2)

        point.fragmentDetails1.setPoint(index1 - dupCount1 - lossCount1)
        point.fragmentDetails2.setPoint(index2 - dupCount2 - lossCount2)

    return points, lostPoints

######################################################
# CountNumLossesAndDuplications
# Parameters:
# Description: Counts the number of duplications and losses before the current index
######################################################
def CountNumLossesAndDuplications(index, losses, dups, strain):
    dupCount = 0
    lossCount = 0
    
    #Number of losses
    if len(losses) > 0:
        for x in  range(0, len(losses)):
            if losses[x].fragmentDetails1.fragmentIndex < index and strain.name == losses[x].genome1Name:
                lossCount += 1
    #Number of duplications
    if len(dups) > 0:
        for x in range(0, len(dups)):
            if dups[x].fragmentDetails1.fragmentIndex < index and strain.name == dups[x].genome1Name:
                dupCount += 1
    return lossCount, dupCount

######################################################
# createBarGraph
# Parameters:
# Description:
######################################################
def createBarGraph(dictionary, title):
    if dictionary != None and len(dictionary) > 0:
        keys = list(dictionary.keys())
        keys.sort()
        
        if len(keys) < 8:
            maxValue = keys[len(keys)-1]
            maxValue += 1
            if (maxValue < 6):
                maxValue = 6
                
            y_pos = []
            for x in range(1, maxValue):
                y_pos.append(x)
            
            performance = []
            for x in range(1, maxValue):
                if x in keys:
                    performance.append(dictionary[x])
                else:
                    performance.append(0)
        else:
            y_pos = np.arange(len(keys))
            performance = []
            for key in keys:
                performance.append(dictionary[key])
        
        #Used for the y-ticks
        val = max(performance)
        ticks = []
        ticks.append(2)
        index = 4        
        while index < val + 2:
            ticks.append(index)
            index += 2
            
        f = plt.figure()
        plt.bar(y_pos, performance, align='center', alpha=0.5)
        if len(keys) < 8:
            plt.xticks(y_pos, y_pos)
        else :
            plt.xticks(y_pos, keys)
        plt.yticks(ticks)
        plt.ylabel('Number of Occurrences')
        plt.xlabel('Size of Occurrence')
        plt.title(title)
        plt.show()
        
        f.savefig("%s.pdf" %(title), bbox_inches='tight')