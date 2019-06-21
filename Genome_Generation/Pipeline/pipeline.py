from genomeGenerator import generateTests
from CompareResultsService import readFiles
import matplotlib.pyplot as plt
import os
import sys
import shutil
import datetime

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
    
    baseCommand = 'python main.py '
    
    currDir = os.getcwd()
    dirList = os.listdir(currDir)
    testingFile = open("tests.txt", "r")
    tests = testingFile.readlines()
    setNumber = 0
    
    totalEventsAveragesList = []
    totalAccuracyAveragesList = []
    strictEventAccuracy = 0.0
    relaxedEventAccuracy = 0.0
    
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
            testSetDir = datetime.datetime.now().strftime("%m-%d-%Y_%H_%M_%S")
            generateTests(testSetDir, tree, maxLength, numOperons, numEvents, probDup, dup_pValue, probLoss, loss_pValue, probInv, inv_pValue, probSub, probTrans, trans_pValue)
            setNumber += 1
#            analyzeTree(tree, testSetDir)
            appCommand = baseCommand + tree + ' ' + testSetDir + ' > appTestingOutput.txt'
            print appCommand
            os.system(appCommand)
            totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents = readFiles(testSetDir)
            
            print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
            if totalEventsExpected > 0:
                strictEventAccuracy = float(totalEventsFound)/float(totalEventsExpected) * 100.0
            if totalGenesExpected > 0:
                relaxedEventAccuracy = float(totalGenesFound)/float(totalGenesExpected) * 100.0
            numEventsAveragesList.append(totalAppEvents)
            accuracyAveragesList.append(strictEventAccuracy)
            
#            count = 0
#            print dirList
#            for testFile in dirList:
#                if testFile.endswith('.pdf'):
#                    print count
#                    count+= 1
#                    shutil.move(os.path.join(currDir, testFile), os.path.join(testSetDir, testFile))
#                elif testFile.startswith('Ancestor'):
#                    shutil.move(os.path.join(currDir, testFile), os.path.join(testSetDir, testFile))
#                elif testFile.startswith('NC_0000'):
#                    shutil.move(os.path.join(currDir, testFile), os.path.join(testSetDir, testFile))
#                elif testFile == "ApplicationOutput.txt":
#                    shutil.move(os.path.join(currDir, testFile), os.path.join(testSetDir, testFile))
#                elif testFile == "generatorOutput.txt":
#                    shutil.move(os.path.join(currDir, testFile), os.path.join(testSetDir, testFile))
                
        totalEventsAveragesList.append(numEventsAveragesList)
        totalAccuracyAveragesList.append(accuracyAveragesList)
        print testSetDir
        
    graphData(totalEventsAveragesList, "Event Averages")
    graphData(totalAccuracyAveragesList, "Accuracy Averages")

def graphData(totalAverages, title):
    averages = []
    xAxis = [25, 50, 75]
    
    for averagesList in totalAverages:
        currentSum = 0.0    
        for average in averagesList:
            currentSum += average
        
        average = currentSum / len(averagesList)
        averages.append(average)
        
    print averages
        
    plt.title(title)
    plt.plot(xAxis, averages, 'o-')
    plt.show()
    

main()