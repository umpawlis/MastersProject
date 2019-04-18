import multiset

####################################
##Sequence Service Functions########
####################################

######################################################
# computeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: Computes the gene content difference between the operons (multiset) ie whats the difference between them while not removing duplicates
######################################################
def computeOperonDifferences(operon1, operon2):
    set1 = multiset.Multiset();
    set2 = multiset.Multiset();
    
    for op in operon1:
        set1.add(op.split('_')[0].strip())
    for op in operon2:
        set2.add(op.split('_')[0].strip())
        
    set3 = set1.symmetric_difference(set2)
    
    return len(set3)