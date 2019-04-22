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

######################################################
# addDuplicationEventsToStrain
# Parameters: stain, duplication sizes, details about the description of the event
# Description: Increments the duplication counter and adds the details to the strain
######################################################
def addDuplicationEventsToStrain(strain, duplicationSizes, duplicationDescription):
    
    if duplicationDescription != None and len(duplicationSizes) > 0:
        strain.duplicationDetails += duplicationDescription
        for x in range(0, len(duplicationSizes)):
            if duplicationSizes[x] in strain.duplicationCounts:
                strain.duplicationCounts[duplicationSizes[x]] += 1
            else:
                strain.duplicationCounts[duplicationSizes[x]] = 1
    return strain

######################################################
# addDeletionEventsToStrain
# Parameters: stain, deletion sizes, details about the description of the event
# Description: Increments the deletion counter and adds the details to the strain
######################################################
def addDeletionEventsToStrain(strain, deletionSizes, deletionDescription):
    
    if deletionDescription != None and len(deletionSizes) > 0:
        strain.deletionDetails += deletionDescription
        for x in range(0, len(deletionSizes)):
            if deletionSizes[x] in strain.deletionCounts:
                strain.deletionCounts[deletionSizes[x]] += 1
            else:
                strain.deletionCounts[deletionSizes[x]] = 1
    return strain