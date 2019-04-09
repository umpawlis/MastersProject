outputFile1 = 'ApplicationOutput.txt'
outputFile2 = 'ApplicationOutput.txt'

######################################################
# parseIntoDictionary
# Parameters:
# Description: parses line into a dictionary
######################################################
def parseIntoDictionary(data):
    dict = {}
    for entry in data:
        entry = entry.split(':')
        size = entry[0]
        occurrences = entry[1]
        dict[size] = occurrences
    return dict

######## Main ########
print('Opening %s %s...' % (outputFile1, outputFile2))
file1 = open(outputFile1, "r")
file2 = open(outputFile2, "r")

if file1.mode == "r" and file2.mode == "r":
    line1 = file1.readline() #Strain
    line2 = file2.readline() #Strain
    
    while line1 and line2:
        #Each loop processes data for one strain
        if line1.strip() == line2.strip():
            currStrain = line1.strip()
            
            line1 = file1.readline() #Deletions
            line2 = file2.readline() #Deletions
            
            if 'Deletions' in line1 and 'Deletions' in line2:
                print('Processing deletions for the following strain: %s' % (currStrain))  
                
                #Parse data
                line1 = line1.replace('Deletions;', '').strip().split(',')
                line2 = line2.replace('Deletions;', '').strip().split(',')
                dict1 = parseIntoDictionary(line1)
                dict2 = parseIntoDictionary(line2)
                
                #TODO calculation
                
                line1 = file1.readline() #Duplications
                line2 = file2.readline() #Duplications
            else:
                print('Error! Expected deletion segments! %s, %s' % (line1, line2))
            
            if 'Duplications' in line1 and 'Duplications' in line2:    
                print('Processing duplications for the following strain: %s' % (currStrain))
                
                #Parse data
                line1 = line1.replace('Duplications;', '').strip().split(',')
                line2 = line2.replace('Duplications;', '').strip().split(',')
                dict1 = parseIntoDictionary(line1)
                dict2 = parseIntoDictionary(line2)
                
                #TODO calculation
                
            else:
                print('Error! Expected duplication segments! %s, %s' % (line1, line2))
        else:
            print('Error! Strains being compared do not match!')
        
        line1 = file1.readline() #Next strain
        line2 = file2.readline() #Next Strain
        
print('Closing files...')
file1.close()
file2.close()
print('Successfully closed files.\nEnd of script...')