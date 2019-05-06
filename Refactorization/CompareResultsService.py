#These indicate the names of the output files
outputFile1 = 'ApplicationOutput.txt'
outputFile2 = 'ApplicationOutput.txt'

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

        line1 = file1.readline() #Strain 1
        line2 = file2.readline() #Strain 2

        while line1 and line2 and len(line1) > 0 and len(line2) > 0:
            if 'Strain' in line1 and 'Strain' in line2: #This is the strain identifier
                strain1 = line1.replace('Strain:', '')
                strain2 = line2.replace('Strain:', '')

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
                        #TODO comparison
                    else:
                        print('Error! This line should be the codon mismatch')
                        return False

                    line1 = file1.readline() #Substitution
                    line2 = file2.readline() #Substitution
                    if 'Substitution' in line1 and 'Substitution' in line2:
                        print('Comparing the substitutions between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the substitutions')
                        return False

                    line1 = file1.readline() #Duplication
                    line2 = file2.readline() #Duplication
                    if 'Duplication' in line1 and 'Duplication' in line2:
                        print('Comparing the duplications between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the duplications')
                        return False

                    line1 = file1.readline() #Deletion
                    line2 = file2.readline() #Deletion
                    if 'Deletion' in line1 and 'Deletion' in line2:
                        print('Comparing the deletions between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the deletions')
                        return False

                    line1 = file1.readline() #Inversion
                    line2 = file2.readline() #Inversion
                    if 'Inversion' in line1 and 'Inversion' in line2:
                        print('Comparing the inversions between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the inversions')
                        return False

                    line1 = file1.readline() #Transposition
                    line2 = file2.readline() #Transposition
                    if 'Transposition' in line1 and 'Transposition' in line2:
                        print('Comparing the transpositions between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the transpositions')
                        return False

                    line1 = file1.readline() #Inverted Transposition
                    line2 = file2.readline() #Inverted Transposition
                    if 'Inverted Transposition' in line1 and 'Inverted Transposition' in line2:
                        print('Comparing the inverted transposition between the strains!')
                        #TODO comparison
                    else:
                        print('Error! This line should be the inverted transpositions')
                        return False
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