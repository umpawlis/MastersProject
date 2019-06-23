import os

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
    
    #Running Duploss
    command = "python " + DUPLOSS_PATH + DUPLOSS_EXEC + " -eiqdt " + inferred + " " + real + " > " + duplossOutFile
    os.system(command)
    
    Tp, Fp, Fn = getTpFpFn(duplossOutFile)
    print "True positive numbers"
    print Tp
    print Fp
    print Fn
    
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
        
    return float(Tp-2), float(Fp), float(Fn)  #removing 2 to TP because we don't want to count the aligned [o] and [t]
    
    
    
##
## Main
##
if __name__ == '__main__':
        
#    real = "< o >, Arg_AGA, [16S, 23S, 5S], Arg_CGU, Ala_GCG, [Leu_CUG, Lys_AAA], Phe_UUC, < t >, -[Met_AUG, Ser_UCG, His_CAU], -Gly_GGA, -Thr_ACG, -[Ile_AUC, Arg_CGG], -[Thr_ACC, Ile_AUC], -[Arg_AGA, Thr_ACC], -[Leu_CUU, Leu_CUA]"
#    inf = "[o],Arg_AGA,Arg_AGA,Leu_UUG,Lys_AAG,16S,23S,5S,Arg_AGG,Arg_CGU,Ala_GCG,Leu_CUG,Ala_GCG,Lys_AAA,Phe_UUC,[t],-Met_AUG,-Ser_UCG,-His_CAU,-Ile_AUC,-Arg_CGG,-Thr_ACC,-Ile_AUC,-Arg_AGA,-Thr_ACC,-Leu_CUU,-Leu_CUA"
    
    real = "< o >, Glu_GAA, Asp_GAC, [16S, 23S, 5S], [Gln_CAA, Phe_UUC], [Gln_CAG, Leu_UUG], [Ser_UCU, Asp_GAU], [Ser_AGC, Ser_UCG, Ser_UCU], [Trp_UGG, Arg_CGA], [Arg_CGC, 16S, Glu_GAA, 23S, 5S], Arg_CGA, [Pro_CCG, Gly_GGA], [Leu_UUA, Arg_AGA], Gly_GGG, [Pro_CCU, Pro_CCU, Ser_UCA, Phe_UUU], [Ala_GCG, Thr_ACG], Pro_CCU, [Val_GUU, Thr_ACG], Ser_UCU, [16S, 23S, 5S], [Pro_CCC, Pro_CCU], [Thr_ACC, His_CAU], [Lys_AAG, Arg_CGC], Arg_CGC, Ala_GCG, [16S, 23S, 5S], Leu_CUC, < t >, -Val_GUU, -[Asp_GAU, Gly_GGA], -[Tyr_UAC, Gln_CAA], -HIs_CAC, -[Trp_UGG, Pro_CCC, Cys_UGU, Ser_AGU], -Arg_CGA, -Gln_CAA, -[Ala_GCC, Gln_CAG], -Phe_UUC, -Leu_CUU, -[Glu_GAA, Asp_GAC], -[Ile_AUU, Ser_UCU, Ser_UCC, Gly_GGU], -Arg_AGG, -[5S, 23S, 16S], -[Pro_CCG, Arg_CGA], -Gln_CAA, -Gly_GGG, -[Gly_GGU, Leu_CUG], -[HIs_CAC, Leu_CUA, 5S, 23S, 16S], -[5S, 23S, 16S], -Thr_ACG, -[Gly_GGC, Ser_AGC, 5S, 23S, 16S], -[Phe_UUU, Gly_GGA], -Gln_CAG"
    inf = "< o >, Glu_GAA, Asp_GAC, [16S, 23S, 5S], [Gln_CAA, Phe_UUC], [Gln_CAG, Leu_UUG], [Ser_UCU, Asp_GAU], [Ser_AGC, Ser_UCG, Ser_UCU], [Trp_UGG, Arg_CGA], [Arg_CGC, 16S, Glu_GAA, 23S, 5S], Arg_CGA, [Pro_CCG, Gly_GGA], [Leu_UUA, Arg_AGA], Gly_GGG, [Pro_CCU, Pro_CCU, Ser_UCA, Phe_UUU], [Ala_GCG, Thr_ACG], Pro_CCU, [Val_GUU, Thr_ACG], Ser_UCU, [16S, 23S, 5S], [Pro_CCC, Pro_CCU], [Thr_ACC, His_CAU], [Lys_AAG, Arg_CGC], Arg_CGC, Ala_GCG, [16S, 23S, 5S], Leu_CUC, < t >, -Val_GUU, -[Asp_GAU, Gly_GGA], -[Tyr_UAC, Gln_CAA], -HIs_CAC, -[Trp_UGG, Pro_CCC, Cys_UGU, Ser_AGU], -Arg_CGA, -Gln_CAA, -[Ala_GCC, Gln_CAG], -Phe_UUC, -Leu_CUU, -[Glu_GAA, Asp_GAC], -[Ile_AUU, Ser_UCU, Ser_UCC, Gly_GGU], -Arg_AGG, -[5S, 23S, 16S], -[Pro_CCG, Arg_CGA], -Gln_CAA, -Gly_GGG, -[Gly_GGU, Leu_CUG], -[HIs_CAC, Leu_CUA, 5S, 23S, 16S], -[5S, 23S, 16S], -Thr_ACG, -[Gly_GGC, Ser_AGC, 5S, 23S, 16S], -[Phe_UUU, Gly_GGA], -Gln_CAG"
    
    recall, precision, fMeasure = compareAnc(inf, real)
    
    print "recall = " + str(recall)
    print "precision = " + str(precision)
    print "f-measure = " + str(fMeasure)