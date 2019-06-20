from genomeGenerator import generateTests
from main import analyzeTree
from CompareResultsService import readFiles
import matplotlib.pyplot as plt
import os
import sys

probDup = 0.0
dup_pValue = 0.0
probLoss = 0.0
loss_pValue = 0.0
probInv = 0.0
inv_pValue = 0.0
probSub = 0.0
probTrans = 0.0
trans_pValue = 0.0

def main():
    global probDup
    global dup_pValue
    global probLoss
    global loss_pValue
    global probInv
    global inv_pValue
    global probSub
    global probTrans
    global trans_pValue
    
    testingFile = open("tests.txt", "r")
    tests = testingFile.readlines()
    setNumber = 0
    
    totalEventsAveragesList = []
    totalAccuracyAveragesList = []
    
    for test in tests:
        numEventsAveragesList = []
        accuracyAveragesList = []
        args = test.split()
        count = 4
        
        tree = args[0]
        maxLength = int(args[1])
        numOperons = int(args[2])
        numEvents = int(args[3])
        
        while count < len(args):
            if args[count] == "-d":
                probDup = float(args[count+1])
                dup_pValue = float(args[count+2])
            elif args[count] == "-l":
                probLoss = float(args[count+1])
                loss_pValue = float(args[count+2])
            elif args[count] == "-i":
                probInv = float(args[count+1])
                inv_pValue = float(args[count+2])
            elif args[count] == "-t":
                probTrans = float(args[count+1])
                trans_pValue = float(args[count+2])
            elif args[count] == "-s":
                probSub = float(args[count+1])
                count -= 1
            count += 3
            
        if probDup + probLoss + probInv + probSub + probTrans != 1.0:
            print "WARNING: Total probability for all events does not equal 1.0. Please change probabilities. Exiting..."
            sys.exit(0)
            
        for i in range(3):
            testSetDir = "Test_Set" + str(setNumber)
            generateTests(testSetDir, tree, maxLength, numOperons, numEvents, probDup, dup_pValue, probLoss, loss_pValue, probInv, inv_pValue, probSub, probTrans, trans_pValue)
            setNumber += 1
            analyzeTree(tree)
            totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents = readFiles()
            
            print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
            strictEventAccuracy = float(totalEventsFound)/float(totalEventsExpected) * 100.0
            relaxedEventAccuracy = float(totalGenesFound)/float(totalGenesExpected) * 100.0
            numEventsAveragesList.append(totalAppEvents)
            accuracyAveragesList.append(strictEventAccuracy)
            
        totalEventsAveragesList.append(numEventsAveragesList)
        totalAccuracyAveragesList.append(accuracyAveragesList)
        
    graphData(totalEventsAveragesList)
    graphData(totalAccuracyAveragesList)

def graphData(totalAverages):
    averages = []
    xAxis = [25, 50]
    
    for averagesList in totalAverages:
        currentSum = 0.0    
        for average in averagesList:
            currentSum += average
        
        average = currentSum / len(averagesList)
        averages.append(average)
        
    plt.title("Data Mapping")
    plt.plot(xAxis, averages, '.-')
    plt.show()
    

main()