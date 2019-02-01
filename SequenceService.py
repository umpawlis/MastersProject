import multiset

####################################
##Sequence Service Functions########
####################################

######################################################
# formatAndComputeOperonDifferences
# Parameters: operon1, operon2 - the two operons to compare
# Description: Creates an array of operons 1 and 2, computes how many genes are unique to 1 operon (no duplicates) and computes the gene content difference between the operons (multiset) ie whats the difference between them while not removing duplicates
######################################################
def formatAndComputeOperonDifferences(operon1, operon2):
    noDuplicatesSet1 = set()
    noDuplicatesSet2 = set()

    set1 = multiset.Multiset();
    set2 = multiset.Multiset();

    operon1 = operon1.replace('-', '')
    operon1 = operon1.replace('[', '')
    operon1 = operon1.replace(']', '')

    operon2 = operon2.replace('-', '')
    operon2 = operon2.replace('[', '')
    operon2 = operon2.replace(']', '')

    operon1List = operon1.split(',')
    operon2List = operon2.split(',')

    noWhiteSpaceOperon1List = []
    noWhiteSpaceOperon2List = []

    for op in operon1List:
        set1.add(op.split('_')[0].strip())
        noDuplicatesSet1.add(op.split('_')[0].strip())
        noWhiteSpaceOperon1List.append(op.strip())

    for op in operon2List:
        set2.add(op.split('_')[0].strip())
        noDuplicatesSet2.add(op.split('_')[0].strip())
        noWhiteSpaceOperon2List.append(op.strip())

    set3 = set1.symmetric_difference(set2)
    noDuplicatesSet3 = noDuplicatesSet1.symmetric_difference(noDuplicatesSet2)

    return len(set3), noWhiteSpaceOperon1List, noWhiteSpaceOperon2List, len(noDuplicatesSet3)

######################################################
# reverseSequence
# Parameters: operon - operon to check
# Description: checks if the operon needs to be reversed
######################################################
def reverseSequence(operon):
    if '-' in operon:
        return True
    return False

######################################################
# formatAllOperons
# Parameters: sequence
# Description: Formats all operons to correct orientation.
######################################################
def formatAllOperons(sequence):
    sequenceList = []
    operonIndexConversions = []
    operonIndex = 0

    for operon1 in sequence:
        reverseOp = reverseSequence(operon1)

        operon1 = operon1.replace('-', '')
        operon1 = operon1.replace('[', '')
        operon1 = operon1.replace(']', '')

        operon1List = operon1.split(',')

        if len(operon1List) > 1:
            noWhiteSpaceOperon1List = []

            for op in operon1List:
                noWhiteSpaceOperon1List.append(op.strip())

            if reverseOp:
                noWhiteSpaceOperon1List.reverse()
            sequenceList.append(noWhiteSpaceOperon1List)
            operonIndexConversions.append(operonIndex)
            operonIndex += 1
        else:
            operonIndexConversions.append(-1)

    return sequenceList, operonIndexConversions