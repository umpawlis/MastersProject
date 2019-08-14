import os
import sys
import subprocess

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
        sys.exit("Error opening file "+filename+': '+str(e))
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
    while genome.count("[t]") > 1: #testing for multiple termini (bug of OrthoAlign)
		print "- one of the genomes have multiple termini, removing the first occurrence"
		genome = genome.replace("[t],", "", 1)
	
    if "< o >" not in genome:  #only does the work if it's in the operon format
        return genome
        
    genome = genome.replace(" ", "") #remove whitespaces
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
    
#    tryAgain = True
#    while tryAgain and count < 2:
    #Running Duploss
    command = "python " + DUPLOSS_PATH + DUPLOSS_EXEC + " -eiqdt " + inferred + " " + real + " > " + duplossOutFile
    os.system(command)
    
    Tp, Fp, Fn, tryAgain = getTpFpFn(duplossOutFile)
#    print "True positive numbers"
#    print Tp
#    print Fp
#    print Fn
    count += 1
    
    if tryAgain:
#        outputFile = open(duplossOutFile, "w+")
        p = subprocess.Popen(['python', DUPLOSS_PATH + DUPLOSS_EXEC, '-eiqdt', inferred, real], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        with open(duplossOutFile, "w+") as f:
            f.write(out)
            f.write(err)
#        outputFile.close()
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
    
    real = "[o],Pro_CCU,Pro_CCU,Ala_GCU,Ala_GCA,Gly_GGG,Ser_AGC,Cys_UGC,Arg_CGA,Ile_AUU,Arg_AGA,Lys_AAA,Pro_CCA,Arg_AGG,Arg_AGG,16S,23S,5S,Gly_GGU,Thr_ACA,Ser_AGU,Pro_CCG,Gly_GGC,16S,23S,5S,16S,23S,5S,Leu_CUC,Met_AUG,Ser_UCC,Ile_AUU,Leu_CUC,Gln_CAG,Ser_AGC,Phe_UUU,Ser_AGC,Arg_AGG,HIs_CAC,Phe_UUC,Asp_GAU,Gln_CAG,Val_GUG,Leu_CUA,Tyr_UAU,-Cys_UGC,-Thr_ACA,-Thr_ACG,-Val_GUA,-5S,-23S,-16S,-Cys_UGU,-HIs_CAC,-Val_GUG,-5S,-23S,-16S,-5S,-23S,-16S,-Val_GUA,-Val_GUU,-Gln_CAG,-Thr_ACA,-Gly_GGU,-Lys_AAG,-Ile_AUU,-Phe_UUU,-5S,-23S,-16S,-Tyr_UAC,-HIs_CAC,-5S,-23S,-16S,-Asp_GAU,-Leu_CUG,-Ile_AUU,-Ser_UCU,-Ser_AGC"
    inf = "[o],Pro_CCU,Pro_CCU,Ala_GCU,Ala_GCA,Gly_GGG,Ser_AGC,Cys_UGC,Arg_CGA,Ile_AUU,Arg_AGA,Lys_AAA,Pro_CCA,Arg_AGG,Arg_AGG,16S,23S,5S,Gly_GGU,Thr_ACA,Ser_AGU,Pro_CCG,Gly_GGC,16S,23S,5S,16S,23S,5S,Leu_CUC,Met_AUG,Ser_UCC,Ile_AUU,Leu_CUC,Gln_CAG,Ser_AGC,Phe_UUU,Ser_AGC,Arg_AGG,HIs_CAC,Phe_UUC,Asp_GAU,Gln_CAG,Val_GUG,Leu_CUA,Tyr_UAU,[t],-Phe_UUU,-Cys_UGC,-Thr_ACA,-Thr_ACG,-Val_GUA,-5S,-23S,-16S,-Cys_UGU,-HIs_CAC,-Val_GUG,-5S,-23S,-16S,-5S,-23S,-16S,-Val_GUA,-Val_GUU,-Gln_CAG,-Thr_ACA,-Lys_AAG,-Ile_AUU,-Phe_UUU,-5S,-23S,-16S,-Tyr_UAC,-HIs_CAC,-5S,-23S,-16S,-Asp_GAU,-Leu_CUG,-Ile_AUU,-Ser_UCU,-Ser_AGC"
        
    recall, precision, fMeasure = compareAnc(inf, real, "")
    
#    print "recall = " + str(recall)
#    print "precision = " + str(precision)
#    print "f-measure = " + str(fMeasure)