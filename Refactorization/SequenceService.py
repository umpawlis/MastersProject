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
                red_x_coord.append(events[i].fragmentDetails1.fragmentIndex)
                red_y_coord.append(events[i].fragmentDetails2.fragmentIndex)
            elif events[i].score == 0:
                green_x_coord.append(events[i].fragmentDetails1.fragmentIndex)
                green_y_coord.append(events[i].fragmentDetails2.fragmentIndex)
            elif events[i].score == 1 or events[i].score == 2:
                yellow_x_coord.append(events[i].fragmentDetails1.fragmentIndex)
                yellow_y_coord.append(events[i].fragmentDetails2.fragmentIndex)
            else:
                orange_x_coord.append(events[i].fragmentDetails1.fragmentIndex)
                orange_y_coord.append(events[i].fragmentDetails2.fragmentIndex)

            #Get all coordinates into a single array
            x_coord.append(events[i].fragmentDetails1.fragmentIndex)
            y_coord.append(events[i].fragmentDetails2.fragmentIndex)
            #print('x-axis: %s, y-axis: %s' %(trackingEvents[i].getGenome1OperonIndex(), trackingEvents[i].getGenome2OperonIndex()))

    #If we have any coordinates to plot, display them
    if len(green_x_coord) > 0 or len(yellow_x_coord) > 0 or len(orange_x_coord) > 0 or len(red_x_coord) > 0:
        f = plt.figure()
        plt.title("Orthologous Operon Mapping")
        plt.plot(green_x_coord, green_y_coord, 'o', color = 'green')
        plt.plot( yellow_x_coord, yellow_y_coord, 'o', color = 'gold')
        plt.plot(orange_x_coord, orange_y_coord, 'o', color = 'orange')
        plt.plot(red_x_coord, red_y_coord, 'o', color = 'red')
        plt.axis([0, len(events)+5, 0, len(events)+5])
        plt.ylabel('Operon Position in %s' % (strain1.name))
        plt.xlabel('Operon Position in %s' % (strain2.name))
        plt.show()
        f.savefig("%s %s.pdf" %(strain1.name, strain2.name), bbox_inches='tight')
    else:
        print('No plot to display!')
    print("x" * 70)
    
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