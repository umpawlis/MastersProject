from genomeGenerator import generateTests
from CompareResultsService import readFiles
from CompareAncestors import compareAnc
from shutil import copy
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

testFolder = ""

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
    global testFolder
    
    cherryTree = False
    if len(sys.argv) < 3:
        print "WARNING: Must provide a file for testing. Exiting..."
        sys.exit(0)
    
    if len(sys.argv) == 5:
        if sys.argv[2] == "-c":
            cherryTree = True
        numRounds = int(sys.argv[3])
        testFolder = sys.argv[4] + "/"
    elif len(sys.argv) == 4:
        if sys.argv[2] == "-c":
            cherryTree = True
        numRounds = int(sys.argv[3])
    elif len(sys.argv) == 3:
        numRounds = int(sys.argv[2])
        
    testFile = sys.argv[1]
    
    
    xAxis = []
    baseCommand = 'python Run_2-SPP_OrthoAlign.py -f '
    
    testingFile = open(testFile, "r")
    testDiff = testingFile.readline().strip()
    tests = testingFile.readlines()
    testingFile.close()
    
    totalEventsAppAveragesList = []
    totalEventsGenAveragesList = []
    totalEventsOrthoAveragesList = []
    totalEventsDupAveragesList = []
    totalStrictAccuracyAveragesList = []
    totalRelaxedAccuracyAveragesList = []
    strictEventAccuracy = 0.0
    relaxedEventAccuracy = 0.0
    
    totalAppFMeasureList = []
    totalOrthoFMeasureList = []
    totalDupFMeasureList = []
    
    for test in tests:
        numEventsAppAveragesList = []
        numEventsGenAveragesList = []
        numEventsOrthoAveragesList = []
        numEventsDupAveragesList = []
        
        strictAccuracyAveragesList = []
        relaxedAccuracyAveragesList = []
        
        appFMeasureList = []
        orthoFMeasureList = []
        dupFMeasureList = []
        
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
            testSetDir = testFolder + datetime.datetime.now().strftime("%m-%d-%Y_%H_%M_%S")
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
            
            if cherryTree:
                appCommand = baseCommand + testSetDir + '/NC_000001/sequence.txt ' + testSetDir + '/NC_000002/sequence.txt ' + testSetDir
                os.system(appCommand)
                
                duplossOutFile = testSetDir + "/duploss.out"
                orthoAlignOutFile = testSetDir + "/orthoAlign.out"
                appRootFile = testSetDir + "/appRoot.txt"
                genRootFile = testSetDir + "/root.txt"
                
                with open(orthoAlignOutFile, "r") as f:
                    line = f.readline()
                    while line:
                        splitted = line.split("ost = ")
                        if len(splitted) > 1:
                            orthoCost = float(splitted[1])
                            line = f.readline()
                        if line.strip() == ">Ancestor:":
                            orthoAncestor = f.readline().strip()
                            break
                        line = f.readline()
                        
                with open(duplossOutFile, "r") as f:
                    line = f.readline()
                    while line:
                        splitted = line.split("ost = ")
                        if len(splitted) > 1:
                            dupCost = float(splitted[1])
                            line = f.readline()
                        if line.strip() == ">Ancestor":
                            dupAncestor = f.readline().strip()
                            break
                        line = f.readline()
                        
                with open(appRootFile, "r") as f:
                    appAncestor = f.readline()
                        
                with open(genRootFile, "r") as f:
                    genAncestor = f.readline()
                
                numEventsOrthoAveragesList.append(orthoCost)
                numEventsDupAveragesList.append(dupCost)
                print orthoCost
                print dupCost
                print orthoAncestor
                print genAncestor
                
                orthoRecall, orthoPrecision, orthofMeasure = compareAnc(orthoAncestor, genAncestor, testSetDir + "/ortho-")
                appRecall, appPrecision, appfMeasure = compareAnc(appAncestor, genAncestor, testSetDir + "/app-")
                dupRecall, dupPrecision, dupfMeasure = compareAnc(dupAncestor, genAncestor, testSetDir + "/dup-")
                appFMeasureList.append(appfMeasure)
                orthoFMeasureList.append(orthofMeasure)
                dupFMeasureList.append(dupfMeasure)
                
                
        totalEventsAppAveragesList.append(numEventsAppAveragesList)
        totalEventsGenAveragesList.append(numEventsGenAveragesList)
        totalEventsOrthoAveragesList.append(numEventsOrthoAveragesList)
        totalEventsDupAveragesList.append(numEventsDupAveragesList)
        
        totalStrictAccuracyAveragesList.append(strictAccuracyAveragesList)
        totalRelaxedAccuracyAveragesList.append(relaxedAccuracyAveragesList)
        
        totalAppFMeasureList.append(appFMeasureList)
        totalOrthoFMeasureList.append(orthoFMeasureList)
        totalDupFMeasureList.append(dupFMeasureList)
        
    outputData(totalEventsAppAveragesList, testFolder + "appEventsData.txt")
    outputData(totalEventsGenAveragesList, testFolder + "genEventsData.txt")
    outputData(totalEventsOrthoAveragesList, testFolder + "orthoEventsData.txt")
    outputData(totalEventsDupAveragesList, testFolder + "dupEventsData.txt")
    
    outputData(totalStrictAccuracyAveragesList, testFolder + "strictAccuracyData.txt")
    outputData(totalRelaxedAccuracyAveragesList, testFolder + "relaxedAccuracyData.txt")
    
    outputData(totalAppFMeasureList, testFolder + "appFMeasureData.txt")
    outputData(totalOrthoFMeasureList, testFolder + "orthoFMeasureData.txt")
    outputData(totalDupFMeasureList, testFolder + "dupFMeasureData.txt")
    
    
    graphData("sAccuracy", totalStrictAccuracyAveragesList, xAxisTitle, xAxis)
    graphData("rAccuracy", totalRelaxedAccuracyAveragesList, xAxisTitle, xAxis)
    if cherryTree:
        graphData("fMeasure", totalAppFMeasureList, xAxisTitle, xAxis, totalAverages3 = totalOrthoFMeasureList, totalAverages4 = totalDupFMeasureList)
        graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList, totalEventsOrthoAveragesList, totalEventsDupAveragesList)
    else:
        graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList)
        
    if testFolder:
        copy(testFile, testFolder)

def graphData(graphType, totalAverages, xAxisTitle, xAxis, totalAverages2 = None, totalAverages3 = None, totalAverages4 = None):
    if graphType == "Events":
        title = "Average Number of Events"
        yAxisTitle = "Number of Events"
    elif graphType == "sAccuracy":
        title = "Average Strict Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracy":
        title = "Average Relaxed Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "fMeasure":
        title = "Average F-measure"
        yAxisTitle = "F-measure"
        
    labels = []
        
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
    labels.append(line1)
    
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
        labels.append(line2)
        
    if totalAverages3 is not None:
        averages3 = []
        print totalAverages3
        for averagesList in totalAverages3:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages3.append(average)
        line3, = plt.plot(xAxis, averages3, 'gP-', label='OrthoAlign')
        labels.append(line3)
        
        
    if totalAverages4 is not None:
        averages4 = []
        print totalAverages4
        for averagesList in totalAverages4:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages4.append(average)
        line4, = plt.plot(xAxis, averages4, 'mX-', label='DupLoss')
        labels.append(line4)
        
    plt.legend(handles=labels)
        
    print averages
    plt.show()
    f.savefig(testFolder + graphType + ".pdf", bbox_inches='tight')
    
def outputData(totalAverages, fileName):
    with open(fileName, 'w+') as f:
        for averagesList in totalAverages:
            for average in averagesList:
                f.write("%s " % (average))
            f.write("\n")

main()