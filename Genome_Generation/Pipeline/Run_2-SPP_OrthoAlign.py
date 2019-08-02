import os
import sys
import getopt
import time

printToConsole = False

### CONSTANTS ###
ORTHOALIGN_PATH =  "OrthoAlign/OrthoAlign/"; ##I recommend using an absolute path
ORTHOALIGN_EXEC = "Aligning"
DUPLOSS_PATH = "2-SPP/"  ##I recommend using an absolute path
DUPLOSS_EXEC = "duploss"

### METHODS ###
def getFirstLineFromFile(filename):
    """
    Return the first line in the file.
    """
    try:
        f = open(filename, 'r')
    except IOError as e:
        sys.exit("Error opening file "+filename+': '+e)
    line = f.readline()
    f.close()
    return line

def remBracketsDistributeNegSign(genome):
    """
    Removes all square brackets (representing operons) and distributes the negative signs in front of operons to all genes inside it.
    Changes <o> and <t> for [o] and [t].
    Returns the new string representing the updated genome.
    """
    negOpSplit = genome.split(",-[")
    
    for i in range(1, len(negOpSplit)):  #skipping the first index --> genomes don't start with -[
        op = negOpSplit[i][:negOpSplit[i].index(']')]  #slice of the operon
        op = op.replace(",", ",-")
        op = '-' + op
        singletons = negOpSplit[i][negOpSplit[i].index(']') + 1:]
        #print "op" + str(i) + " = " + op
        negOpSplit[i] = op + singletons
        
    genome = ",".join(negOpSplit)
    genome = genome.replace("[", "")
    genome = genome.replace("]", "")
    genome = genome.replace("<o>", "[o]")
    return genome.replace("<t>", "[t]")
    
def getCostAndAncestorFromOutFile(outfile):
    """
    Opens either a duploss or orthoAlign output file, finds the cost and returns it.
    """
    try:
        f = open(outfile, 'r')
    except IOError as e:
        sys.exit("Error opening file "+outfile+': '+e)
        
    line = f.readline()
    while line:
        splitted = line.split("ost = ")
        if len(splitted) > 1:
            cost = float(splitted[1])
            break
        line = f.readline()
    
    line = f.readline()
    while line:
        if ">Ancestor" in line:  #the next line should be the ancestor
            ancestor = f.readline()
            break
        line = f.readline()
    
    f.close()
    return cost, ancestor
    

##
## Main
##
if __name__ == '__main__':

    usage = "usage: Run_2-SPP_OrthoAlign -hf genome1 genome2 destinationFolderName\n"+\
            " -f genome1 and genome2 are filenames instead of genome strings.\n"\
            " -h print this message.\n"\

    try:
        opts, args = getopt.getopt(sys.argv[1:], "fh")
    except getopt.GetoptError:
        sys.exit(usage)

    input_files = False
    
    for (opt,val) in opts:
        if(opt == "-h"):
            sys.exit(usage)
        if(opt == "-f"):
            input_files = True
    
    if len(args) != 3:
        sys.exit(usage)
        
        
    if(input_files):
        file1 = args[0]
        file2 = args[1]
        if(not os.path.exists(file1)):
            sys.exit("Error: "+file1+" does not exist.\n")
        if(not os.path.exists(file2)):
            sys.exit("Error: "+file2+" does not exist.\n")
        genome1 = getFirstLineFromFile(file1)
        genome2 = getFirstLineFromFile(file2)
    else:
        genome1 = args[0]
        genome2 = args[1]
        
    #Output filenames with output prefix
    duplossOutFile = args[2] + "/duploss.out"
    orthoAlignOutFile = args[2] + "/orthoAlign.out"
        
    #Removing all whitespaces
    genome1 = genome1.replace(" ", "")
    genome2 = genome2.replace(" ", "")
    
    #Removing square brackets and distributing negative signs inside operons, updating origin and terminus
    genome1 = remBracketsDistributeNegSign(genome1)
    genome2 = remBracketsDistributeNegSign(genome2)
        
    #print "genome 1 = " + genome1 + "\ngenome 2 = " + genome2
    
    #Running OrthoAlign
    command = "java -classpath " + ORTHOALIGN_PATH + " " + ORTHOALIGN_EXEC + " -dt " + genome1 + " " + genome2 + " > " + orthoAlignOutFile
    orthoAlignStartTime = time.time()
    os.system(command)
    orthoAlignRunTime = time.time() - orthoAlignStartTime
    
    #Running Duploss
    command = "python " + DUPLOSS_PATH + DUPLOSS_EXEC + " -eiqdt " + genome1 + " " + genome2 + " > " + duplossOutFile
    duplossStartTime = time.time()
    os.system(command)
    duplossRunTime = time.time() - duplossStartTime
    
    orthoCost, orthoAncestor = getCostAndAncestorFromOutFile(orthoAlignOutFile)
    duplossCost, duplossAncestor = getCostAndAncestorFromOutFile(duplossOutFile)
    
    if printToConsole:
        print "Duploss cost = " + str(duplossCost)
        print "Duploss ancestor =    " + duplossAncestor
        print "OrthoAlign cost = " + str(orthoCost)
        print "OrthoAlign ancestor = " + orthoAncestor
        
    fileDirectory = args[2].split("/")
    runtimePath = "/".join(fileDirectory[:-1])
    with open(runtimePath + "/OrthoRuntimes.txt", "a+") as runtimeFile:
        runtimeFile.write("%f " % (orthoAlignRunTime))
        
    with open(runtimePath + "/DuplossRuntimes.txt", "a+") as runtimeFile:
        runtimeFile.write("%f " % (duplossRunTime))