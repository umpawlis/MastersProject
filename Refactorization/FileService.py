####################################################
######File Service Functions Refactorization########
####################################################

######################################################
# createFile
# Parameters:
# Description: Creates a file where the output will be stored
######################################################
def createFile(fileName):
    print('Creating file %s...' % (fileName));
    file = open(fileName, "w+")
    file.close()