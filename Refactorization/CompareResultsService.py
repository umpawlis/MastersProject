#These indicate the names of the output files
outputFile1 = 'ApplicationOutput.txt'
outputFile2 = 'ApplicationOutput.txt'

######################################################
# inversionTranspositionComparison
# Parameters:
# Description: Compares results of inversions, transpositions, and inverted transpositions
######################################################
def inversionTranspositionComparison(data1, data2):
    percentage = 0
    dict1 = {}
    dict2 = {}

    #Parse the data
    regions1 = data1.split('|') #An inversion/transposition fragment ie entire piece that was transposed
    regions2 = data2.split('|') #An inversion/transposition fragment ie entire piece that was transposed
    for region in regions1:
        operons = region.split(';')
        for operon in operons:
            genes = operon.split(',')
            for gene in genes:
                data = gene.split(' ')
                if len(data) == 2:
                    dict1[data[1]] = data[0]
    for region in regions2:
        operons = region.split(';')
        for operon in operons:
            genes = operon.split(',')
            for gene in genes:
                data = gene.split(' ')
                if len(data) == 2:
                    dict2[data[1]] = data[0]
    #Compute a percentage
    keys = dict1.keys()
    count = 0
    for key in keys:
        if key in dict2 and dict2[key] == dict1[key]: #A correctly identified event
            count += 1

    if count == 0 and len(dict2) == 0:
        return 100
    else:
        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return percentage

######################################################
# duplicationDeletionComparison
# Parameters:
# Description: Compares results of codon mismatches and substitutions
######################################################
def duplicationDeletionComparison(data1, data2):
    percentage = 0

    #Parse the data
    segments1 = data1.split(';')
    segments2 = data2.split(';')
    dict1 = {}
    dict2 = {}
    for segment in segments1:
        genes = segment.split(',')
        for gene in genes:
            data = gene.split(' ')
            if len(data) == 2:
                dict1[data[1]] = data[0]
    for segment in segments2:
        genes = segment.split(',')
        for gene in genes:
            data = gene.split(' ')
            if len(data) == 2:
                dict2[data[1]] = data[0]

    #Compute a percentage
    keys = dict1.keys()
    count = 0
    for key in keys:
        if key in dict2 and dict2[key] == dict1[key]: #A correctly identified event
            count += 1

    if count == 0 and len(dict2) == 0:
        return 100
    else:
        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return percentage

######################################################
# codonMismatchSubstitutionComparison
# Parameters:
# Description: Compares results of codon mismatches and substitutions
######################################################
def codonMismatchSubstitutionComparison(data1, data2):
    percentage = 0
    dict1 = {}
    dict2 = {}

    #Parse the data
    array1 = data1.split(';')
    array2 = data2.split(';')
    for entry in array1:
        data = entry.split(' ')
        if len(data) == 2:
            dict1[data[1]] = data[0]
    for entry in array2:
        data = entry.split(' ')
        if len(data) == 2:
            dict2[data[1]] = data[0]

    #Compute a percentage
    keys = dict1.keys()
    count = 0
    for key in keys:
        if key in dict2 and dict2[key] == dict1[key]: #A correctly identified event
            count += 1

    if count == 0 and len(dict2) == 0:
        return 100
    else:
        percentage = (count/len(dict2)) * 100 #Number of correct events divided by the total events from simulator

    return percentage

######################################################
# readFiles
# Parameters:
# Description: Reads two files and compares two files
######################################################
def readFiles():
    newickTree1 = ''
    newickTree2 = ''
    print('Opening %s %s...' % (outputFile1, outputFile2))
    file1 = open(outputFile1, "r")
    file2 = open(outputFile2, "r")

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
                    print('Comparing the following strain between the two files: %s' % (strain1))

                    line1 = file1.readline() #Codon mismatch
                    line2 = file2.readline() #Codon mismatch
                    if 'Codon Mismatch' in line1 and 'Codon Mismatch' in line2:
                        print('Comparing the codon mismatches between the strains!')
                        line1 = line1.replace('Codon Mismatch:', '')
                        line2 = line2.replace('Codon Mismatch:', '')
                        result = codonMismatchSubstitutionComparison(line1, line2)
                        print('The result of the Codon Mismatches is: %s percent' % (result))
                    else:
                        print('Error! This line should be the codon mismatch')
                        return False

                    line1 = file1.readline() #Substitution
                    line2 = file2.readline() #Substitution
                    if 'Substitution' in line1 and 'Substitution' in line2:
                        print('Comparing the substitutions between the strains!')
                        line1 = line1.replace('Substitution:', '')
                        line2 = line2.replace('Substitution:', '')
                        result = codonMismatchSubstitutionComparison(line1, line2)
                        print('The result of the Substitutions is: %s percent' % (result))
                    else:
                        print('Error! This line should be the substitutions')
                        return False

                    line1 = file1.readline() #Duplication
                    line2 = file2.readline() #Duplication
                    if 'Duplication' in line1 and 'Duplication' in line2:
                        print('Comparing the duplications between the strains!')
                        line1 = line1.replace('Duplication:', '')
                        line2 = line2.replace('Duplication:', '')
                        result = duplicationDeletionComparison(line1, line2)
                        print('The result of the Duplications is: %s percent' % (result))
                    else:
                        print('Error! This line should be the duplications')
                        return False

                    line1 = file1.readline() #Deletion
                    line2 = file2.readline() #Deletion
                    if 'Deletion' in line1 and 'Deletion' in line2:
                        print('Comparing the deletions between the strains!')
                        line1 = line1.replace('Deletion:', '')
                        line2 = line2.replace('Deletion:', '')
                        result = duplicationDeletionComparison(line1, line2)
                        print('The result of the Deletions is: %s percent' % (result))
                    else:
                        print('Error! This line should be the deletions')
                        return False

                    line1 = file1.readline() #Inversion
                    line2 = file2.readline() #Inversion
                    if 'Inversion' in line1 and 'Inversion' in line2:
                        print('Comparing the inversions between the strains!')
                        line1 = line1.replace('Inversion:', '')
                        line2 = line2.replace('Inversion:', '')
                        result = inversionTranspositionComparison(line1, line2)
                        print('The result of the Inversion is: %s percent' % (result))
                    else:
                        print('Error! This line should be the inversions')
                        return False

                    line1 = file1.readline() #Transposition
                    line2 = file2.readline() #Transposition
                    if 'Transposition' in line1 and 'Transposition' in line2:
                        print('Comparing the transpositions between the strains!')
                        line1 = line1.replace('Transposition:', '')
                        line2 = line2.replace('Transposition:', '')
                        result = inversionTranspositionComparison(line1, line2)
                        print('The result of the Transposition is: %s percent' % (result))
                    else:
                        print('Error! This line should be the transpositions')
                        return False

                    line1 = file1.readline() #Inverted Transposition
                    line2 = file2.readline() #Inverted Transposition
                    if 'Inverted Transposition' in line1 and 'Inverted Transposition' in line2:
                        print('Comparing the inverted transposition between the strains!')
                        line1 = line1.replace('Inverted Transposition:', '')
                        line2 = line2.replace('Inverted Transposition:', '')
                        result = inversionTranspositionComparison(line1, line2)
                        print('The result of the Inverted Transposition is: %s percent' % (result))
                    else:
                        print('Error! This line should be the inverted transpositions')
                        return False
            print('\n')
            line1 = file1.readline() #Strain 1
            line2 = file2.readline() #Strain 2
    else:
        print('Unable to process output files!')
    print('Closing files...')
    file1.close()
    file2.close()
    print('Successfully closed files.')
    return True

######## Main ########
if readFiles():
    print('Successfully processed the output files')
else:
    print('Error! An error has occured while processing the files!')
print('End of script...')