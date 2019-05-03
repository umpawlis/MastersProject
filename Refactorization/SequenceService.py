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
def createDotPlot(events, strain1, strain2):
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

    print("x" * 70)
    for i in range(0, len(events)):
        if events[i].technique == 'Global Alignment' or events[i].technique == 'Local Alignment':
            #Assign the coords to the appropriate array, the index represents the position of the operon with respect to the genome
            if events[i].technique == 'Local Alignment':
                red_x_coord.append(events[i].fragmentDetails1.point)
                red_y_coord.append(events[i].fragmentDetails2.point)
            elif events[i].score == 0 and (events[i].fragmentDetails1.description == 'Terminus' or events[i].fragmentDetails1.description == 'Origin'):
                black_x_coord.append(events[i].fragmentDetails1.point)
                black_y_coord.append(events[i].fragmentDetails2.point)
            elif events[i].score == 0:
                green_x_coord.append(events[i].fragmentDetails1.point)
                green_y_coord.append(events[i].fragmentDetails2.point)
            elif events[i].score == 1 or events[i].score == 2:
                yellow_x_coord.append(events[i].fragmentDetails1.point)
                yellow_y_coord.append(events[i].fragmentDetails2.point)
            else:
                orange_x_coord.append(events[i].fragmentDetails1.point)
                orange_y_coord.append(events[i].fragmentDetails2.point)

            #Get all coordinates into a single array
            x_coord.append(events[i].fragmentDetails1.point)
            y_coord.append(events[i].fragmentDetails2.point)
            #print('x-axis: %s, y-axis: %s' %(events[i].fragmentDetails1.point, events[i].fragmentDetails2.point))

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
        f.savefig("%s %s.pdf" %(strain1.name, strain2.name), bbox_inches='tight')
    else:
        print('No plot to display!')
    print("x" * 70)


######################################################
# adjustOperonIndexesForPlot
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

        if point.fragmentDetails1.fragmentIndex == point.fragmentDetails2.fragmentIndex: #The operon's position is conserved keep as is
            point.fragmentDetails1.setPoint(point.fragmentDetails1.fragmentIndex)
            point.fragmentDetails2.setPoint(point.fragmentDetails2.fragmentIndex)
        else: #These points are different, need to adjust it according to the losses and duplications
            index1 = point.fragmentDetails1.fragmentIndex
            index2 = point.fragmentDetails2.fragmentIndex

            count1 = CountNumLossesAndDuplications(index1, lostPoints, dup1, strain1)
            count2 = CountNumLossesAndDuplications(index2, lostPoints, dup2, strain2)

            point.fragmentDetails1.setPoint(index1 - count1)
            point.fragmentDetails2.setPoint(index2 - count2)

    return points, lostPoints

######################################################
# CountNumLossesAndDuplications
# Parameters:
# Description: Counts the number of duplications and losses before the current index
######################################################
def CountNumLossesAndDuplications(index, losses, dups, strain):
    count = 0
    #Number of losses
    if len(losses) > 0:
        for x in  range(0, len(losses)):
            if losses[x].fragmentDetails1.fragmentIndex < index and strain.name == losses[x].genome1Name:
                count += 1
    #Number of duplications
    if len(dups) > 0:
        for x in  range(0, len(dups)):
            if dups[x].fragmentDetails1.fragmentIndex < index and strain.name == dups[x].genome1Name:
                count += 1
    return count

######################################################
# createBarGraph
# Parameters:
# Description:
######################################################
def createBarGraph(dictionary, title):
    if dictionary != None and len(dictionary) > 0:
        keys = list(dictionary.keys())
        keys.sort()

        y_pos = np.arange(len(keys))

        performance = []
        for key in keys:
            performance.append(dictionary[key])

        plt.bar(y_pos, performance, align='center', alpha=0.5)
        plt.xticks(y_pos, keys)
        plt.ylabel('Number of Occurrences')
        plt.xlabel('Size of Occurrence')
        plt.title(title)
        plt.show()