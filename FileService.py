
####################################
######File Service Functions########
####################################

######################################################
# createFile
# Parameters:
# Description: Creates a file where the output will be stored
######################################################
def createFile(fileName):
    print('Creating file %s...' % (fileName));
    file = open(fileName, "w+")
    file.close()
    
######################################################
# appendToFile
# Parameters:
# Description: Appends content to file
######################################################
def appendToFile(fileName, content):
    print('Appending to file %s the following content, %s...' % (fileName, content))
    file = open(fileName, "a")
    file.write(content)
    file.close()