from genomeGenerator import generateTests
from CompareResultsService import readFiles
import matplotlib.pyplot as plt
import os
import sys
import subprocess
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
    
    if len(sys.argv) != 3:
        print "WARNING: Must provide a file for testing. Exiting..."
        sys.exit(0)
        
    testFile = sys.argv[1]
    numRounds = int(sys.argv[2])
    
    xAxis = []
    baseCommand = 'python main.py '
    
    testingFile = open(testFile, "r")
    testDiff = testingFile.readline().strip()
    tests = testingFile.readlines()
    testingFile.close()
    
    totalEventsAppAveragesList = []
    totalEventsGenAveragesList = []
    totalStrictAccuracyAveragesList = []
    totalRelaxedAccuracyAveragesList = []
    strictEventAccuracy = 0.0
    relaxedEventAccuracy = 0.0
    
    for test in tests:
        numEventsAppAveragesList = []
        numEventsGenAveragesList = []
        strictAccuracyAveragesList = []
        relaxedAccuracyAveragesList = []
        args = test.strip().split()
        count = 4
        
        tree = args[0]
        maxLength = int(args[1])
        numOperons = int(args[2])
        numEvents = int(args[3])
        
        if testDiff == "Genes":
            xAxisTitle = "Size of Genome"
            xAxis.append(maxLength)
        elif testDiff == "Events":
            xAxisTitle = "Number of Events per Branch"
            xAxis.append(numEvents)
        elif testDiff == "Tree":
            xAxisTitle = "Size of Tree"
            if len(tree) > 4:
                xAxis.append(int(tree[4]))
            else:
                print "WARNING: Tree file must be in format tree#*.dnd where # is the number of leaves. Exiting..."
                sys.exit(0)
        
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
            
        for i in range(numRounds):
            testSetDir = datetime.datetime.now().strftime("%m-%d-%Y_%H_%M_%S")
            generateTests(testSetDir, tree, maxLength, numOperons, numEvents, probDup, dup_pValue, probLoss, loss_pValue, probInv, inv_pValue, probSub, probTrans, trans_pValue)
#            analyzeTree(tree, testSetDir)
#            appCommand = baseCommand + tree + ' ' + testSetDir + ' > ' + testSetDir + '/appTestingOutput.txt'
#            os.system(appCommand)
#            subprocess.Popen(appCommand, shell=True).wait()
            p = subprocess.Popen(['python', 'main.py', tree, testSetDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            with open(testSetDir + '/appTestingOutput.txt', "w+") as f:
                f.write(out)
                f.write(err)
            
            totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents = readFiles(testSetDir)
            
            print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected, totalAppEvents))
            if totalEventsExpected > 0:
                strictEventAccuracy = float(totalEventsFound)/float(totalEventsExpected) * 100.0
            else:
                strictEventAccuracy = 0.0
            if totalGenesExpected > 0:
                relaxedEventAccuracy = float(totalGenesFound)/float(totalGenesExpected) * 100.0
            else:
                relaxedEventAccuracy = 0.0
                
            numEventsAppAveragesList.append(totalAppEvents)
            numEventsGenAveragesList.append(totalEventsExpected)
            strictAccuracyAveragesList.append(strictEventAccuracy)
            relaxedAccuracyAveragesList.append(relaxedEventAccuracy)
                
        totalEventsAppAveragesList.append(numEventsAppAveragesList)
        totalEventsGenAveragesList.append(numEventsGenAveragesList)
        totalStrictAccuracyAveragesList.append(strictAccuracyAveragesList)
        totalRelaxedAccuracyAveragesList.append(relaxedAccuracyAveragesList)
        print testSetDir
        
    outputData(totalEventsAppAveragesList, "appEventsData.txt")
    outputData(totalEventsGenAveragesList, "genEventsData.txt")
    outputData(totalStrictAccuracyAveragesList, "strictAccuracyData.txt")
    outputData(totalRelaxedAccuracyAveragesList, "relaxedAccuracyData.txt")
    graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList)
    graphData("sAccuracy", totalStrictAccuracyAveragesList, xAxisTitle, xAxis)
    graphData("rAccuracy", totalRelaxedAccuracyAveragesList, xAxisTitle, xAxis)

def graphData(graphType, totalAverages, xAxisTitle, xAxis, totalAverages2 = None):
    if graphType == "Events":
        title = "Average Number of Events"
        yAxisTitle = "Number of Events"
    elif graphType == "sAccuracy":
        title = "Average Strict Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracy":
        title = "Average Relaxed Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
    f = plt.figure()
    plt.title(title)
    plt.ylabel(yAxisTitle)
    plt.xlabel(xAxisTitle)
    plt.grid(True)
        
    averages = [] 
    print totalAverages
    for averagesList in totalAverages:
        currentSum = 0.0    
        for average in averagesList:
            currentSum += average
        
        average = currentSum / len(averagesList)
        averages.append(average)
    line1, = plt.plot(xAxis, averages, 'o-', label='Application')
    
    averages2 = []
    if totalAverages2 is not None:
        print totalAverages2
        for averagesList in totalAverages2:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages2.append(average)
        line2, = plt.plot(xAxis, averages2, 'r^-', label='Generator')
        plt.legend(handles=[line1, line2])
    else:
        plt.legend(handles=[line1])
        
    print averages
    plt.show()
    f.savefig(graphType + ".pdf", bbox_inches='tight')
    
def outputData(totalAverages, fileName):
    with open(fileName, 'w+') as f:
        for averagesList in totalAverages:
            for average in averagesList:
                f.write("%s " % (average))
            f.write("\n")

main()