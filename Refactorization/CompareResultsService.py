#These indicate the names of the output files
outputFile1 = 'ApplicationOutput.txt'
outputFile2 = 'ApplicationOutput.txt'

######## Main ########
print('Opening %s %s...' % (outputFile1, outputFile2))
file1 = open(outputFile1, "r")
file2 = open(outputFile2, "r")

newickTree1 = ''
newickTree2 = ''

if file1.mode == "r" and file2.mode == "r":
    newickTree1 = file1.readline() #Newick tree 1
    newickTree2 = file2.readline() #Newick tree 2

else:
    print('Unable to process output files!')
    
print('Closing files...')
file1.close()
file2.close()
print('Successfully closed files.\nEnd of script...')