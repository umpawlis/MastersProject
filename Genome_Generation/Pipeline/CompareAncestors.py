import os
import sys

### CONSTANTS ###
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

def cleanUpGenomes(genome):
    """
    Removes whitespaces.
    Removes all square brackets (representing operons) and distributes the negative signs in front of operons to all genes inside it.
    Changes <o> and <t> for [o] and [t].
    Returns the new string representing the updated genome.
    """
    if "< o >" not in genome:  #only does the work if it's in the operon format
        return genome
        
    genome = genome.replace(" ", "") #remove whitespaces
    negOpSplit = genome.split(",-[")
    
    for i in range(1, len(negOpSplit)):  #skipping the first index --> genomes don't start with -[
        op = negOpSplit[i][:negOpSplit[i].index(']')]  #slice of the operon
        op = op.replace(",", ",-")
        op = '-' + op
        #print "op" + str(i) + " = " + op
        negOpSplit[i] = op
        
    genome = ",".join(negOpSplit)
    genome = genome.replace("[", "")
    genome = genome.replace("]", "")
    genome = genome.replace("<o>", "[o]")
    return genome.replace("<t>", "[t]")
    
    
def compareAnc(inferred, real, outputPath):
    """
    Compare the inferred ancestor (from app, duploss or orthoAlign) with the real (generated) ancestor.
    Returns the recall (sensitivity), precision, and F-measure (harmonic mean of both).
    """
    inferred = cleanUpGenomes(inferred)
    real = cleanUpGenomes(real)
    
    #Output filenames with output prefix
    duplossOutFile = outputPath + "compareAnc-duploss.out"
    count = 0
    
    tryAgain = True
    while tryAgain and count < 2:
        #Running Duploss
        command = "python " + DUPLOSS_PATH + DUPLOSS_EXEC + " -eiqdt " + inferred + " " + real + " > " + duplossOutFile
        os.system(command)
        
        Tp, Fp, Fn, tryAgain = getTpFpFn(duplossOutFile)
        print "True positive numbers"
        print Tp
        print Fp
        print Fn
        count += 1
        print "Count: " + str(count)
    
    if Tp != 0:
        recall = Tp / (Tp + Fn)
        precision = Tp / (Tp + Fp)
        fMeasure = 2 * (recall * precision) / (recall + precision)
    else:
        recall = 0
        precision = 0
        fMeasure = 0
    
    #print "Tp = " + str(Tp)
    #print "Fp = " + str(Fp)
    #print "Fn = " + str(Fn)
    
    return recall, precision, fMeasure
    
    
def getTpFpFn(outfile):
    """
    Opens a duploss output file, in which the X genome is the inferred ancestor and the Y genome is the real ancestor.
    Finds and returns:
    TP = number of matches
    FP = number of extra genes
    FN = number of missing genes
    """
    try:
        f = open(outfile, 'r')
    except IOError as e:
        sys.exit("Error opening file "+outfile+': '+e)
        
    Tp = 0
    Fp = 0
    Fn = 0
    tryAgain = False
    startCounting = False  #start counting TP, FP, FN
    line = f.readline()
    while line:
        if not startCounting and "========" in line:
            startCounting = True
        
        elif startCounting:
            splitted = line.split(" --- ")
            if len(splitted) < 2: #it means we're done parsing the alignment
                break
            
            if "." in splitted[0]: #it means a gene is missing!
                Fn += 1
            elif "." in splitted[1]: #it means there's an extra gene!
                Fp += 1
            else: #it means we have a match!
                Tp += 1
            
        line = f.readline()
    if not startCounting:
        Tp = 2
        print "Duploss file is incomplete for %s. Difference between genomes was too great." % (outfile)
        tryAgain = True
        
    return float(Tp-2), float(Fp), float(Fn), tryAgain   #removing 2 to TP because we don't want to count the aligned [o] and [t]
    
    
    
##
## Main
##
if __name__ == '__main__':
        
#    real = "< o >, Arg_AGA, [16S, 23S, 5S], Arg_CGU, Ala_GCG, [Leu_CUG, Lys_AAA], Phe_UUC, < t >, -[Met_AUG, Ser_UCG, His_CAU], -Gly_GGA, -Thr_ACG, -[Ile_AUC, Arg_CGG], -[Thr_ACC, Ile_AUC], -[Arg_AGA, Thr_ACC], -[Leu_CUU, Leu_CUA]"
#    inf = "[o],Arg_AGA,Arg_AGA,Leu_UUG,Lys_AAG,16S,23S,5S,Arg_AGG,Arg_CGU,Ala_GCG,Leu_CUG,Ala_GCG,Lys_AAA,Phe_UUC,[t],-Met_AUG,-Ser_UCG,-His_CAU,-Ile_AUC,-Arg_CGG,-Thr_ACC,-Ile_AUC,-Arg_AGA,-Thr_ACC,-Leu_CUU,-Leu_CUA"
    
#    real = ""
#    inf = ""
    
    real = "[o],Arg_AGG,Gly_GGA,Gly_GGU,Val_GUA,Ile_AUA,Ala_GCC,Ser_UCG,Ser_UCA,Asp_GAC,Ser_UCU,Ser_UCA,Gln_CAA,Ile_AUA,Ala_GCG,Gly_GGC,Phe_UUC,Tyr_UAU,Gly_GGG,Ala_GCA,Asp_GAU,Ser_UCU,Phe_UUU,Leu_CUU,Ile_AUU,16S,23S,5S,Ser_UCA,16S,23S,5S,Ala_GCA,Gly_GGC,Ala_GCA,Ser_UCU,Val_GUC,Ala_GCG,Leu_UUG,Leu_CUU,Pro_CCG,Ile_AUU,Val_GUG,Gly_GGG,Ala_GCG,Ile_AUA,-Ser_UCU,-5S,-23S,-16S,-Leu_UUG,-Ala_GCA,-Val_GUU,-Arg_AGA,-Ser_AGU,-Cys_UGU,-Gly_GGA,-Arg_CGC,-Ile_AUA,-Ser_UCA,-Arg_CGG,-Tyr_UAC,-Pro_CCU,-Trp_UGG,-Pro_CCC,-Gln_CAA,-Ala_GCC,-Leu_UUG,-Leu_CUC,-5S,-23S,-16S,-Arg_CGC,-His_CAU,-Gly_GGG,-Arg_CGA,-Arg_AGA,-Leu_UUA,-HIs_CAC,-Arg_CGU,-Ala_GCG,-Arg_CGC,-5S,-23S,-16S"
    inf = "[o],Arg_AGG,Gly_GGA,Gly_GGU,Val_GUA,Ile_AUA,Ala_GCC,Ser_UCG,Ser_UCA,Asp_GAC,Ser_UCU,Ser_UCA,Gly_GGA,Ala_GCG,Gly_GGC,Gln_CAA,Ile_AUA,Phe_UUC,Tyr_UAU,Gly_GGG,Ala_GCA,Asp_GAU,Ser_UCU,Leu_CUU,Ile_AUU,16S,23S,5S,Phe_UUU,Ser_UCA,16S,23S,5S,Ala_GCA,Gly_GGC,Ala_GCA,Ser_UCU,Val_GUC,Ala_GCG,Leu_UUG,Leu_CUU,Pro_CCG,Ile_AUU,Val_GUG,Gly_GGG,16S,23S,5S,[t],-Leu_UUG,-Ala_GCA,-Val_GUU,-Arg_AGA,-Ile_AUA,-Ala_GCG,-Ser_AGU,-Cys_UGU,-Gly_GGA,-Arg_CGC,-Ile_AUA,-Ser_UCA,-Arg_CGG,-Tyr_UAC,-Pro_CCU,-Trp_UGG,-Pro_CCC,-Gln_CAA,-Ala_GCC,-Leu_UUG,-Leu_CUC,-5S,-23S,-16S,-Arg_CGC,-His_CAU,-Gly_GGG,-Arg_CGA,-Arg_AGA,-Leu_UUA,-HIs_CAC,-Arg_CGU,-Ala_GCG,-Arg_CGC,-5S,-23S,-16S"
        
    recall, precision, fMeasure = compareAnc(inf, real, "")
    
    print "recall = " + str(recall)
    print "precision = " + str(precision)
    print "f-measure = " + str(fMeasure)