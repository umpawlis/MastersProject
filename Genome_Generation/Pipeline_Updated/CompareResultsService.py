import numpy as np
import matplotlib.pyplot as plt

#These indicate the names of the output files
#outputFile1 = 'ApplicationOutput.txt'
#outputFile2 = 'generatorOutput.txt'

minSize = 100
maxSize = 0
medianSize = 0
sizeSum = 0
totalNumGenes = 0
totalGenes = []

def resetGlobals():
    global minSize
    global maxSize
    global sizeSum
    global medianSize
    global totalNumGenes
    global totalGenes
    
    minSize = 0
    maxSize = 0
    sizeSum = 0
    medianSize = 0
    totalNumGenes = 0
    totalGenes = []

######################################################
# constructDistributionGraph
# Parameters:
# Description: Constructs a distribution graph
######################################################
def constructDistributionGraph(dictionary1, dictionary2, title):
    print('Constructing distribution graph...')
    maxSize1 = 0
    maxSize2 = 0
    
    if dictionary1 != None and len(dictionary1) > 0:
        keys = dictionary1.keys()
        for key in keys:
            if int(key) > maxSize1:
                maxSize1 = int(key)
    if dictionary2 != None and len(dictionary2) > 0:
        keys = dictionary2.keys()
        for key in keys:
            if int(key) > maxSize2:
                maxSize2 = int(key)   
    maxSize = max(maxSize1, maxSize2) + 1 #Make sure to include the last key!
    
    if maxSize != 0:
        x = np.arange(maxSize)
        y1 = []
        y2 = []
        
        for i in range(0, maxSize):
            if str(i) in dictionary1:
                y1.append(int(dictionary1[str(i)]))
            else:
                y1.append(0)
                
            if str(i) in dictionary2:
                y2.append(int(dictionary2[str(i)]))
            else:
                y2.append(0)
        
        plt.bar(x + 0.00, y1, color = 'b', width = 0.25)
        plt.bar(x + 0.25, y2, color = 'g', width = 0.25)
        plt.xticks(x)
        
        plt.ylabel('Number of Events')
        plt.xlabel('Size of Event')
        plt.title(title)
        plt.show()
    
    print('Done constructing graph.')

######################################################
# parseSizeDistribution
# Parameters:
# Description: Parses the size distribution into a dictionary
######################################################
def parseSizeDistribution(line):
    dictionary = {}
    array = line.split(',')
    
    if len(array) > 0:
        for item in array:
            data = item.strip()
            dataArray = data. split()
            
            if len(dataArray) != 4:
#                print('Error! Should be 4 elements for the size distribution')
                return {}
            else:
                size = dataArray[1]
                count = dataArray[3]
                if size in dictionary:
                    dictionary[size] += count
                else:
                    dictionary[size] = count
    return dictionary

######################################################
# inversionTranspositionComparison
# Parameters:
# Description: Compares results of inversions, transpositions, and inverted transpositions
######################################################
def inversionTranspositionComparison(data1, data2, outputFile):
    global minSize
    global maxSize
    global sizeSum
    global totalNumGenes
    global totalGenes
    
    percentage = 0
    dict1 = {}
    dict2 = {}
    negativePositionsList1 = []
    negativePositionsList2 = []

    #Parse the data
    regions1 = data1.strip().split('|') #An inversion/transposition fragment ie entire piece that was transposed
    outputFile.write('|'.join(regions1) + "\n")
    regions2 = data2.strip().split('|') #An inversion/transposition fragment ie entire piece that was transposed
    outputFile.write('|'.join(regions2)  + "\n")
    
    numEventsFound = 0
    numEventsExpected = len(regions2) - 1
    numGenesFound = 0
    numGenesExpected = 0
    numAppEvents = len(regions1) - 1
    
    count = -1
    for region in regions1:
        operons = region.split(';')
        for operon in operons:
            genes = operon.split(', ')
            for gene in genes:
                data = gene.split(' ')
                if len(data) == 2:
                    if data[1] != '-1':
                        dict1[data[1]] = data[0]
                    else:
                        negativePositionsList1.append(data[0])
    for region in regions2:
        eventFound = False
        if region != '':
#            outputFile.write("Reversing Section\n")
#            outputFile.write(region + '\n')
#            reversedRegion = getReversed(region) +';'
#            outputFile.write(reversedRegion + '\n')
#            if region in regions1:
#                numEventsFound += 1
            if flexCompareEvents(region, regions1):
                numEventsFound += 1
                eventFound = True
            operons = region.split(';')
            for operon in operons:
                if operon != '':
                    genes = operon.split(', ')
                    numGenesExpected += len(genes)
                    if eventFound:
                        sizeSum += len(genes)
                        totalGenes.append(len(genes))
                        if len(genes) < minSize:
                            minSize = len(genes)
                        if len(genes) > maxSize:
                            maxSize = len(genes)
                    for gene in genes:
                        data = gene.split(' ')
                        if len(data) == 2:
                            if data[1] != '-1':
                                dict2[data[1]] = data[0]
                            else:
                                dict2[str(count)] = 'XXX'
                                count -= 1
                                negativePositionsList2.append(data[0])
                                
#    print dict1
#    print dict2
#    print negativePositionsList1
#    print negativePositionsList2
    
    #Compute a percentage
    keys = dict2.keys()
    count = 0
    for key in keys:
        if int(key) < 0:
            if count < len(negativePositionsList2):
                if negativePositionsList2[count] in negativePositionsList1:
                    numGenesFound += 1
            count += 1
        else:
            for newKey in range(int(key)-2, int(key)+3):
                if str(newKey) in dict1 and dict2[key] == dict1[str(newKey)]: #A correctly identified event
                    numGenesFound += 1
                    break

#    if count == 0 and len(dict2) == 0:
#        return 100
#    else:
#        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return (numEventsFound, numEventsExpected, numGenesFound, numGenesExpected, numAppEvents)

######################################################
# duplicationDeletionComparison
# Parameters:
# Description: Compares results of codon mismatches and substitutions
######################################################
def duplicationDeletionComparison(data1, data2, outputFile):
    global minSize
    global maxSize
    global sizeSum
    global totalNumGenes
    global totalGenes
    
    percentage = 0

    #Parse the data
    segments1 = data1.strip().split(';')
    outputFile.write('|'.join(segments1) + "\n")
    segments2 = data2.strip().split(';')
    outputFile.write('|'.join(segments2) + "\n")
    dict1 = {}
    dict2 = {}
    negativePositionsList1 = []
    negativePositionsList2 = []
    
    numEventsFound = 0
    numEventsExpected = len(segments2) - 1
    numGenesFound = 0
    numGenesExpected = 0
    numAppEvents = len(segments1) - 1
    
    count = -1
    for segment in segments1:
        genes = segment.split(', ')
        for gene in genes:
            data = gene.split(' ')
            if len(data) == 2:
                if data[1] != '-1':
                    dict1[data[1]] = data[0]
                else:
                    negativePositionsList1.append(data[0])
    for segment in segments2:
        eventFound = False
        if segment != '':
#            if ',' in segment:
#                outputFile.write("Reversing Section\n")
#                outputFile.write(segment + '\n')
#                reversedSegment = getReversed(segment)
#                outputFile.write(reversedSegment + '\n')
#                if reversedSegment in segments1:
#                    numEventsFound += 1
#            if segment in segments1:
#                numEventsFound += 1
            if flexCompareEvents(segment, segments1):
                numEventsFound += 1
                eventFound = True
            genes = segment.split(', ')
            numGenesExpected += len(genes)
            if eventFound:
                sizeSum += len(genes)
                totalGenes.append(len(genes))
                if len(genes) < minSize:
                    minSize = len(genes)
                if len(genes) > maxSize:
                    maxSize = len(genes)
            for gene in genes:
                data = gene.split(' ')
                if len(data) == 2:
                    if data[1] != '-1':
                        dict2[data[1]] = data[0]
                    else:
                        dict2[str(count)] = 'XXX'
                        count -= 1
                        negativePositionsList2.append(data[0])
                        
#    print dict1
#    print dict2
#    print negativePositionsList1
#    print negativePositionsList2

    #Compute a percentage
    keys = dict2.keys()
    count = 0
    for key in keys:
        if int(key) < 0:
            if count < len(negativePositionsList2):
                if negativePositionsList2[count] in negativePositionsList1:
                    numGenesFound += 1
            count += 1
        else:
            for newKey in range(int(key)-2, int(key)+3):
                if str(newKey) in dict1 and dict2[key] == dict1[str(newKey)]: #A correctly identified event
                    numGenesFound += 1
                    break

#    if count == 0 and len(dict2) == 0:
#        return 100
#    else:
#        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return (numEventsFound, numEventsExpected, numGenesFound, numGenesExpected, numAppEvents)

######################################################
# codonMismatchSubstitutionComparison
# Parameters:
# Description: Compares results of codon mismatches and substitutions
######################################################
def codonMismatchSubstitutionComparison(data1, data2, outputFile):
    global minSize
    global maxSize
    global sizeSum
    global totalNumGenes
    global totalGenes
    
    percentage = 0
    dict1 = {}
    dict2 = {}
    negativePositionsList1 = []
    negativePositionsList2 = []

    #Parse the data
    array1 = data1.strip().split(';')
    outputFile.write('|'.join(array1) + "\n")
    array2 = data2.strip().split(';')
    outputFile.write('|'.join(array2) + "\n")
    
    numEventsFound = 0
    numEventsExpected = len(array2) - 1
    numAppEvents = len(array1) - 1
    
    count = -1
    for entry in array1:
        data = entry.split(' ')
        if len(data) == 2:
            if data[1] != '-1':
                dict1[data[1]] = data[0]
            else:
                negativePositionsList1.append(data[0])
    for entry in array2:
        if entry != '':
            data = entry.split(' ')
            if len(data) == 2:
                if data[1] != '-1':
                    dict2[data[1]] = data[0]
                else:
                    dict2[str(count)] = 'XXX'
                    count -= 1
                    negativePositionsList2.append(data[0])
                    
#    print dict1
#    print dict2
#    print negativePositionsList1
#    print negativePositionsList2

    #Compute a percentage
    keys = dict2.keys()
    count = 0
    for key in keys:
        if int(key) < 0:
            if count < len(negativePositionsList2):
                if negativePositionsList2[count] in negativePositionsList1:
                    numEventsFound += 1
                    sizeSum += 1
                    totalGenes.append(1)
                    if minSize > 1:
                        minSize = 1
                    if maxSize < 1:
                        maxSize = 1
            count += 1
        else:
            for newKey in range(int(key)-2, int(key)+3):
                if str(newKey) in dict1 and dict2[key] == dict1[str(newKey)]: #A correctly identified event
                    numEventsFound += 1
                    sizeSum += 1
                    totalGenes.append(1)
                    if minSize > 1:
                        minSize = 1
                    if maxSize < 1:
                        maxSize = 1
                    break

#    if count == 0 and len(dict2) == 0:
#        return 100
#    else:
#        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return (numEventsFound, numEventsExpected, numAppEvents)

def flexCompareEvents(event, compareList):
    formattedEvent = []
    if ';' in event:
        operons = event.split(';')
        for operon in operons:
            if operon != '':
                genes = operon.split(', ')
                formattedEvent += genes
    else:
        genes = event.split(', ')
        formattedEvent += genes
#    print formattedEvent
    
    equal = False
    for compareEvent in compareList:
        if compareEvent != '':
            formattedCompareEvent = []
            if ';' in compareEvent:
                operons = compareEvent.split(';')
                for operon in operons:
                    if operon != '':
                        genes = operon.split(', ')
                        formattedCompareEvent += genes
            else:
                genes = compareEvent.split(', ')
                formattedCompareEvent += genes
                
#            print formattedCompareEvent
                
            if len(formattedEvent) == len(formattedCompareEvent):
                for i in range(len(formattedEvent)):
                    if '< t >' in formattedEvent[i] or '< o >' in formattedEvent[i]:
                        gene1 = formattedEvent[i][:5]
                        position1 = int(formattedEvent[i][6:])
                    else:
                        data1 = formattedEvent[i].split(' ')
                        gene1 = data1[0]
                        position1 = int(data1[1])
                        
                    if '< t >' in formattedCompareEvent[i] or '< o >' in formattedCompareEvent[i]:
                        gene2 = formattedCompareEvent[i][:5]
                        position2 = int(formattedCompareEvent[i][6:])
                    else:
                        data2 = formattedCompareEvent[i].split(' ')
                        gene2 = data2[0]
                        position2 = int(data2[1])
                    
                    if gene1 == gene2:
#                        print gene1 + " " + gene2
#                        print str(position1)
#                        print range(position2-2, position2+3)
                        if position1 in range(position2-2, position2+3):
                            if i == len(formattedEvent)-1:
                                equal = True
                                return equal
                        else:
                            break
                    else:
                        break
    return equal
                

def getReversed(section):
    operons = section.strip().split(';')
    operons.reverse()
    
    reversedSection = []
    for operon in operons:
        if operon != "":
            if ', ' in operon:
                genes = operon.strip().split(', ')
                genes.reverse()
                reversedSection.append(', '.join(genes))
            else:
                reversedSection.append(operon)
                
    return ';'.join(reversedSection)

######################################################
# getTotalEvents
# Parameters:
# Description: computes the total number of deletions or duplications
######################################################
def getTotalEvents(line):
    result = 0
    array = line.split(',')
    
    if len(array) > 0:
        for data in array:
            processedLine = (data.replace('size:', '')).replace('count:', '')
            numbers = processedLine.split()
            if len(numbers) == 2:
                size = int(numbers[0])
                count = int(numbers[1])
                result += (size * count)
            else:
                print('Error! There should be two numbers!')
    return result

######################################################
# readFiles
# Parameters:
# Description: Reads two files and compares two files
######################################################
def readFiles(fileDir, outputFile1, outputFile2, prefix):
    global medianSize
    newickTree1 = ''
    newickTree2 = ''
    outputFile = open(fileDir+ "/" + prefix + "comparisonOutput.txt", "w+")
    outputFile.write('Opening %s %s...\n' % (outputFile1, outputFile2))
    file1 = open(fileDir+ "/" + outputFile1, "r")
    file2 = open(fileDir+ "/" + outputFile2, "r")
    
    totalEventsFound = 0
    totalEventsExpected = 0
    totalGenesFound = 0
    totalGenesExpected = 0
    totalAppEvents = 0
    
    resetGlobals()
    
    #Track total counts for each event type
    duplicationTotals = [0, 0, 0, 0]
    lossTotals = [0, 0, 0, 0]
    inversionTotals = [0, 0, 0, 0]
    transpositionTotals = [0, 0, 0, 0]
    substitutionTotals = [0, 0, 0, 0]
    
    if file1.mode == "r" and file2.mode == "r":
        newickTree1 = file1.readline() #Newick tree 1
        newickTree2 = file2.readline() #Newick tree 2
        
        if newickTree1 == newickTree2:
            print('Trees match!')

        line1 = file1.readline() #Strain 1
        line2 = file2.readline() #Strain 2

        while line1 and line2 and len(line1) > 0 and len(line2) > 0:
            if 'Strain' in line1 and 'Strain' in line2: #This is the strain identifier
                strain1 = line1.replace('Strain:', '').strip()
                strain2 = line2.replace('Strain:', '').strip()

                #Indicate whether this is a matching strain
                if strain1 != strain2:
                    print('Error! The strain names do not match!')
                    return False
                else:
                    outputFile.write('Comparing the following strain between the two files: %s\n' % (strain1))

                    line1 = file1.readline() #Codon mismatch
                    line2 = file2.readline() #Codon mismatch
                    line2 = file2.readline() #Substitution instead because generator does not produce codon mismatches
                    if 'Codon Mismatch' in line1 and 'Substitution' in line2:
                        outputFile.write('Comparing the codon mismatches between the strains!\n')
                        codonLine = line1.replace('Codon Mismatch:', '')
#                        line2 = line2.replace('Codon Mismatch:', '')
                        line2 = line2.replace('Substitution:', '')
#                        result = codonMismatchSubstitutionComparison(line1, line2, outputFile)
#                        
#                        totalEventsFound += result[0]
#                        substitutionTotals[0] += result[0]
#                        
#                        totalEventsExpected += result[1]
#                        substitutionTotals[1] += result[1]
#                        
#                        totalGenesFound += result[0]
#                        substitutionTotals[2] += result[0]
#                        
#                        totalGenesExpected += result[1]
#                        substitutionTotals[3] += result[1]
#                        
#                        totalAppEvents += result[2]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Codon Mismatches is: %s percent' % (result))
                    else:
                        print('Error! This line should be the codon mismatch')
                        return False

                    line1 = file1.readline() #Substitution
#                    line2 = file2.readline() #Substitution
#                    if 'Substitution' in line1 and 'Substitution' in line2:
                    if 'Substitution' in line1:
                        outputFile.write('Comparing the substitutions between the strains!\n')
                        line1 = line1.replace('Substitution:', '').strip()
                        line1 = line1 + codonLine
#                        line2 = line2.replace('Substitution:', '')
                        result = codonMismatchSubstitutionComparison(line1, line2, outputFile)
                        
                        totalEventsFound += result[0]
                        substitutionTotals[0] += result[0]
                        
                        totalEventsExpected += result[1]
                        substitutionTotals[1] += result[1]
                        
                        totalGenesFound += result[0]
                        substitutionTotals[2] += result[0]
                        
                        totalGenesExpected += result[1]
                        substitutionTotals[3] += result[1]
                        
                        totalAppEvents += result[2]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Substitutions is: %s percent' % (result))
                    else:
                        print('Error! This line should be the substitutions')
                        return False

                    line1 = file1.readline() #Duplication
                    line2 = file2.readline() #Duplication
                    if 'Duplication' in line1 and 'Duplication' in line2:
                        outputFile.write('Comparing the duplications between the strains!\n')
                        line1 = line1.replace('Duplication:', '')
                        line2 = line2.replace('Duplication:', '')
                        result = duplicationDeletionComparison(line1, line2, outputFile)
                        
                        totalEventsFound += result[0]
                        duplicationTotals[0] += result[0]
                        
                        totalEventsExpected += result[1]
                        duplicationTotals[1] += result[1]
                        
                        totalGenesFound += result[2]
                        duplicationTotals[2] += result[2]
                        
                        totalGenesExpected += result[3]
                        duplicationTotals[3] += result[3]
                        
                        totalAppEvents += result[4]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Duplications is: %s percent' % (result))
                    else:
                        print('Error! This line should be the duplications')
                        return False

                    line1 = file1.readline() #Deletion
                    line2 = file2.readline() #Deletion
                    if 'Deletion' in line1 and 'Deletion' in line2:
                        outputFile.write('Comparing the deletions between the strains!\n')
                        line1 = line1.replace('Deletion:', '')
                        line2 = line2.replace('Deletion:', '')
                        result = duplicationDeletionComparison(line1, line2, outputFile)
                        
                        totalEventsFound += result[0]
                        lossTotals[0] += result[0]
                        
                        totalEventsExpected += result[1]
                        lossTotals[1] += result[1]
                        
                        totalGenesFound += result[2]
                        lossTotals[2] += result[2]
                        
                        totalGenesExpected += result[3]
                        lossTotals[3] += result[3]
                        
                        totalAppEvents += result[4]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Deletions is: %s percent' % (result))
                    else:
                        print('Error! This line should be the deletions')
                        return False

                    line1 = file1.readline() #Inversion
                    line2 = file2.readline() #Inversion
                    if 'Inversion' in line1 and 'Inversion' in line2:
                        outputFile.write('Comparing the inversions between the strains!\n')
                        line1 = line1.replace('Inversion:', '')
                        line2 = line2.replace('Inversion:', '')
                        result = inversionTranspositionComparison(line1, line2, outputFile)
                        
                        totalEventsFound += result[0]
                        inversionTotals[0] += result[0]
                        
                        totalEventsExpected += result[1]
                        inversionTotals[1] += result[1]
                        
                        totalGenesFound += result[2]
                        inversionTotals[2] += result[2]
                        
                        totalGenesExpected += result[3]
                        inversionTotals[3] += result[3]
                        
                        totalAppEvents += result[4]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Inversion is: %s percent' % (result))
                    else:
                        print('Error! This line should be the inversions')
                        return False

                    line1 = file1.readline() #Transposition
                    line2 = file2.readline() #Transposition
                    if 'Transposition' in line1 and 'Transposition' in line2:
                        outputFile.write('Comparing the transpositions between the strains!\n')
                        line1 = line1.replace('Transposition:', '')
                        line2 = line2.replace('Transposition:', '')
                        result = inversionTranspositionComparison(line1, line2, outputFile)
                        
                        totalEventsFound += result[0]
                        transpositionTotals[0] += result[0]
                        
                        totalEventsExpected += result[1]
                        transpositionTotals[1] += result[1]
                        
                        totalGenesFound += result[2]
                        transpositionTotals[2] += result[2]
                        
                        totalGenesExpected += result[3]
                        transpositionTotals[3] += result[3]
                        
                        totalAppEvents += result[4]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Transposition is: %s percent' % (result))
                    else:
                        print('Error! This line should be the transpositions')
                        return False
                    
                    line1 = file1.readline() #Inverted Transposition
                    line2 = file2.readline() #Inverted Transposition
                    if 'Inverted Transposition' in line1 and 'Inverted Transposition' in line2:
                        outputFile.write('Comparing the inverted transposition between the strains!\n')
                        line1 = line1.replace('Inverted Transposition:', '')
                        line2 = line2.replace('Inverted Transposition:', '')
                        result = inversionTranspositionComparison(line1, line2, outputFile)                        
                        
                        totalEventsFound += result[0]
                        transpositionTotals[0] += result[0]                        
                        
                        totalEventsExpected += result[1]
                        transpositionTotals[1] += result[1]                        
                        
                        totalGenesFound += result[2]
                        transpositionTotals[2] += result[2]                        
                        
                        totalGenesExpected += result[3]
                        transpositionTotals[3] += result[3]
                        
                        totalAppEvents += result[4]
                        outputFile.write('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s\n' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
#                        print('The result of the Inverted Transposition is: %s percent' % (result))
                    else:
                        print('Error! This line should be the inverted transpositions')
                        return False
                    
            elif 'Total Deletions' in line1 and 'Total Deletions' in line2:
                outputFile.write('Comparing total deletions between files!\n')
                line1 = line1.replace('Total Deletions:', '').strip()
                line2 = line2.replace('Total Deletions:', '').strip()
                count1 = getTotalEvents(line1)
                count2 = getTotalEvents(line2)
                if count2 > 0:
                    accuracyRate = (count1/count2) * 100
                    outputFile.write('Accuracy rate for deletions was %s %%\n' % (accuracyRate))
                
                line1 = file1.readline() #Total duplications
                line2 = file2.readline() #Total duplications
                if 'Total Duplications' in line1 and 'Total Duplications' in line2:
                    outputFile.write('Comparing total duplications between files!\n')
                    line1 = line1.replace('Total Duplications:', '').strip()
                    line2 = line2.replace('Total Duplications:', '').strip()
                    count1 = getTotalEvents(line1)
                    count2 = getTotalEvents(line2)
                    if count2 > 0:
                        accuracyRate = (count1/count2) * 100
                        outputFile.write('Accuracy rate for deletions was %s %%\n' % (accuracyRate))
                else:
                    print('Error! Expected total duplications!')
                    return False
                
                line1 = file1.readline() #Size distribution of inversions
                line2 = file2.readline() #Size distribution of inversions
                if 'Size Distribution of Inversions:' in line1 and 'Size Distribution of Inversions:' in line2:
#                    print('Comparing size distribution for inversions')
                    line1 = line1.replace('Size Distribution of Inversions:', '').strip()
                    line2 = line2.replace('Size Distribution of Inversions:', '').strip()
                    dict1 = parseSizeDistribution(line1)
                    dict2 = parseSizeDistribution(line2)
#                    constructDistributionGraph(dict1, dict2, 'Size Distribution for Inversions')
                else:
                    print('Error! Expected size distribution of inversions')
                
                line1 = file1.readline() #Size distribution of transpositions
                line2 = file2.readline() #Size distribution of transpositions
                if 'Size Distribution of Transpositions:' in line1 and 'Size Distribution of Transpositions:' in line2:
#                    print('Comparing size distribution for transpositions')
                    line1 = line1.replace('Size Distribution of Transpositions:', '').strip()
                    line2 = line2.replace('Size Distribution of Transpositions:', '').strip()
                    dict1 = parseSizeDistribution(line1)
                    dict2 = parseSizeDistribution(line2)
#                    constructDistributionGraph(dict1, dict2, 'Size Distribution for Transpositions')
                else:
                    print('Error! Expected size distribution of transpositions')
                    
                line1 = file1.readline() #Size distribution of inverted transpositions
                line2 = file2.readline() #Size distribution of inverted transpositions
                if 'Size Distribution of Inverted Transpositions:' in line1 and 'Size Distribution of Inverted Transpositions:' in line2:
#                    print('Comparing size distribution for inverted transpositions')
                    line1 = line1.replace('Size Distribution of Inverted Transpositions:', '').strip()
                    line2 = line2.replace('Size Distribution of Inverted Transpositions:', '').strip()
                    dict1 = parseSizeDistribution(line1)
                    dict2 = parseSizeDistribution(line2)
#                    constructDistributionGraph(dict1, dict2, 'Size Distribution for Inverted Transpositions')
                else:
                    print('Error! Expected size distribution of inverted transpositions')
                    
                line1 = file1.readline() #Total inversions
                line2 = file2.readline() #Total inversions
                if 'Total Inversions' in line1 and 'Total Inversions' in line2:
#                    print('Comparing total inversions between files!')
                    line1 = line1.replace('Total Inversions:', '').strip()
                    line2 = line2.replace('Total Inversions:', '').strip()
                    count1 = int(line1) #Count of events
                    count2 = int(line2) #Count of events
                    if count2 > 0:
                        accuracyRate = (count1/count2) * 100
#                        print('Accuracy rate for inversions was %s %%' % (accuracyRate))
                else:
                    print('Error! Expected total inversions!')
                    return False
                
                line1 = file1.readline() #Total transpositions
                line2 = file2.readline() #Total transpositions
                if 'Total Transpositions' in line1 and 'Total Transpositions' in line2:
#                    print('Comparing total transpositions between files!')
                    line1 = line1.replace('Total Transpositions:', '').strip()
                    line2 = line2.replace('Total Transpositions:', '').strip()
                    count1 = int(line1) #Count of events
                    count2 = int(line2) #Count of events
                    if count2 > 0:
                        accuracyRate = (count1/count2) * 100
#                        print('Accuracy rate for transpositions was %s %%' % (accuracyRate))
                else:
                    print('Error! Expected total transpositions!')
                    return False
                
                line1 = file1.readline() #Total inverted transpositions
                line2 = file2.readline() #Total inverted transpositions
                if 'Total Inverted Transpositions' in line1 and 'Total Inverted Transpositions' in line2:
#                    print('Comparing total inverted transpositions between files!')
                    line1 = line1.replace('Total Inverted Transpositions:', '').strip()
                    line2 = line2.replace('Total Inverted Transpositions:', '').strip()
                    count1 = int(line1) #Count of events
                    count2 = int(line2) #Count of events
                    if count2 > 0:
                        accuracyRate = (count1/count2) * 100
#                        print('Accuracy rate for inverted transpositions was %s %%' % (accuracyRate))
                        
                    #Run to end of file
                    line1 = file1.readline() #Codon mismatch
                    line2 = file2.readline() #Codon mismatch
                    line1 = file1.readline() #subs
                    line2 = file2.readline() #subs
                    line1 = file1.readline() #time
                    line2 = file2.readline() #time
                else:
                    print('Error! Expected total inverted transpositions!')
                    return False
   
            print('\n')
            line1 = file1.readline() #Strain 1 or Totals at the end
            line2 = file2.readline() #Strain 2 or Totals at the end
    else:
        print('Unable to process output files!')
        
    outputFile.write('Duplication Events Found: %s Duplication Events Expected: %s Duplication Genes Found: %s Duplication Genes Expected: %s\n' % (duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3]))
    outputFile.write('Deletion Events Found: %s Deletion Events Expected: %s Deletion Genes Found: %s Deletion Genes Expected: %s\n' % (lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3]))
    outputFile.write('Inversion Events Found: %s Inversion Events Expected: %s Inversion Genes Found: %s Inversion Genes Expected: %s\n' % (inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3]))
    outputFile.write('Transposition Events Found: %s Transposition Events Expected: %s Transposition Genes Found: %s Transposition Genes Expected: %s\n' % (transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3]))
#    outputFile.write('Inverted Transposition Events Found: %s Inverted Transposition Events Expected: %s Inverted TranspositionInverted Transposition Genes Found: %s Inverted Transposition Genes Expected: %s\n' % (invertedTranspositionTotals[0], invertedTranspositionTotals[1], invertedTranspositionTotals[2], invertedTranspositionTotals[3]))
    outputFile.write('Substitution Events Found: %s Substitution Events Expected: %s Substitution Genes Found: %s Substitution Genes Expected: %s\n' % (substitutionTotals[0], substitutionTotals[1], substitutionTotals[2], substitutionTotals[3]))
    
#    print('Closing files...')
    outputFile.close()
    file1.close()
    file2.close()
#    print('Successfully closed files.')
    
    directories = fileDir.split("/")
    dataFileDir = "/".join(directories[:-1])
    
#    avgEventSize = float(sizeSum) / float(totalEventsFound)
    with open(dataFileDir + "/" + prefix + "EventSizeData.txt", "a+") as dataFile:
        dataFile.write("%f " % sizeSum)
    with open(dataFileDir + "/" + prefix + "EventCountData.txt", "a+") as dataFile:
        dataFile.write("%f " % totalEventsFound)
    with open(dataFileDir + "/" + prefix + "EventMinData.txt", "a+") as dataFile:
        if minSize == 100:
            dataFile.write("0 ")
        else:
            dataFile.write("%f " % minSize)
    with open(dataFileDir + "/" + prefix + "EventMaxData.txt", "a+") as dataFile:
        dataFile.write("%f " % maxSize)
    if len(totalGenes) > 0:
        totalGenes.sort()
        if (len(totalGenes) % 2) == 1:
            medianSize = totalGenes[len(totalGenes)/2]
        else:
            medianSize = (totalGenes[len(totalGenes)/2] + totalGenes[(len(totalGenes)/2) - 1]) / 2
    with open(dataFileDir + "/" + prefix + "EventMedianData.txt", "a+") as dataFile:
        dataFile.write("%f " % medianSize)
        
    return totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals, substitutionTotals

######## Main ########
#if readFiles("compareTest", 'ApplicationOutput.txt', 'generatorOutput.txt', 'app-'):
#    print('Successfully processed the output files')
#else:
#    print('Error! An error has occured while processing the files!')
#print('End of script...')