
outputFile1 = 'ApplicationOutput.txt'
outputFile2 = 'ApplicationOutput.txt'

######## Main ########
print('Starting output comparison script...\nThe following files will be compared: %s, %s\nOpening files...' %(outputFile1, outputFile2))

file1 = open(outputFile1, "r")
file2 = open(outputFile2, "r")

if file1.mode == "r" and file2.mode == "r":
    print('Successfully opened the files')
    line1 = file1.readline()
    line2 = file2.readline()
    
    while line1 and line2:
        
        #TODO
        
        line1 = file1.readline()
        line2 = file2.readline()

print('Closing files...')
file1.close()
file2.close()
print('Successfully closed files.\nEnd of script...')