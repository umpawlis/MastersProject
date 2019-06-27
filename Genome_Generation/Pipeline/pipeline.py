from genomeGenerator import generateTests
from CompareResultsService import readFiles
from CompareAncestors import compareAnc
from ListEvents import outputEvents
from shutil import copy
import matplotlib.pyplot as plt
import os
import sys
import subprocess
import datetime

printToConsole = True

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
    neighbour = False
    equalEvents = False
    if len(sys.argv) < 3:
        print "WARNING: Must provide a file for testing. Exiting..."
        sys.exit(0)
    
    if len(sys.argv) == 5:
        if 'c' in sys.argv[2]:
            cherryTree = True
            if 'n' in sys.argv[2]:
                neighbour = True
        if 'e' in sys.argv[2]:
            equalEvents = True
            
        numRounds = int(sys.argv[3])
        testFolder = sys.argv[4] + "/"
    elif len(sys.argv) == 4:
        if 'c' in sys.argv[2]:
            cherryTree = True
            if 'n' in sys.argv[2]:
                neighbour = True
            numRounds = int(sys.argv[3])
        if 'e' in sys.argv[2]:
            equalEvents = True
            numRounds = int(sys.argv[3])
            
        if not cherryTree and not equalEvents:
            numRounds = int(sys.argv[2])
            testFolder = sys.argv[3] + "/"
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
    
    totalEventsAppNeighbourAveragesList = []
    
    totalStrictAppAccuracyAveragesList = []
    totalRelaxedAppAccuracyAveragesList = []
    totalStrictOrthoAccuracyAveragesList = []
    totalRelaxedOrthoAccuracyAveragesList = []
    totalStrictDupAccuracyAveragesList = []
    totalRelaxedDupAccuracyAveragesList = []
    # strictEventAccuracy = 0.0
    # relaxedEventAccuracy = 0.0
    totalStrictAppNeighbourAccuracyAveragesList = []
    totalRelaxedAppNeighbourAccuracyAveragesList = []
    
    totalAppFMeasureList = []
    totalOrthoFMeasureList = []
    totalDupFMeasureList = []
    
    totalAppNeighbourFMeasureList = []
    
    for test in tests:
        numEventsAppAveragesList = []
        numEventsGenAveragesList = []
        numEventsOrthoAveragesList = []
        numEventsDupAveragesList = []
        numEventsAppNeighbourAveragesList = []
        
        strictAppAccuracyAveragesList = []
        relaxedAppAccuracyAveragesList = []
        strictOrthoAccuracyAveragesList = []
        relaxedOrthoAccuracyAveragesList = []
        strictDupAccuracyAveragesList = []
        relaxedDupAccuracyAveragesList = []
        
        strictAppNeighbourAccuracyAveragesList = []
        relaxedAppNeighbourAccuracyAveragesList = []
        
        appFMeasureList = []
        orthoFMeasureList = []
        dupFMeasureList = []
        
        appNeighbourFMeasureList = []
        
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
        
        basePValue = 0.0
        while count < len(args):
            if args[count] == "-d":
                probDup = float(args[count+1])
                dup_pValue = float(args[count+2])
                basePValue = dup_pValue
            elif args[count] == "-l":
                probLoss = float(args[count+1])
                loss_pValue = float(args[count+2])
                basePValue = loss_pValue
            elif args[count] == "-i":
                probInv = float(args[count+1])
                inv_pValue = float(args[count+2])
                basePValue = inv_pValue
            elif args[count] == "-t":
                probTrans = float(args[count+1])
                trans_pValue = float(args[count+2])
                basePValue = trans_pValue
            elif args[count] == "-s":
                probSub = float(args[count+1])
                count -= 1
            count += 3
        
        if testDiff == "pValue":
            if basePValue != 0.0:
                xAxisTitle = "Value of Geometric Sampling Parameter"
                xAxis.append(basePValue)
            else:
                print "WARNING: Must have atleast one pValue for the test. Exiting..."
                sys.exit(0)
            
        if probDup + probLoss + probInv + probSub + probTrans != 1.0:
            print "WARNING: Total probability for all events does not equal 1.0. Please change probabilities. Exiting..."
            sys.exit(0)
            
        for i in range(numRounds):
            testSetDir = testFolder + datetime.datetime.now().strftime("%m-%d-%Y_%H_%M_%S")
            if neighbour:
                tree = 'tree2LeafNeighbour.dnd'
            generateTests(testSetDir, tree, maxLength, numOperons, numEvents, probDup, dup_pValue, probLoss, loss_pValue, probInv, inv_pValue, probSub, probTrans, trans_pValue, equalEvents)
            if neighbour:
                tree = args[0]
#            analyzeTree(tree, testSetDir)
#            appCommand = baseCommand + tree + ' ' + testSetDir + ' > ' + testSetDir + '/appTestingOutput.txt'
#            os.system(appCommand)
#            subprocess.Popen(appCommand, shell=True).wait()
            p = subprocess.Popen(['python', 'main.py', tree, testSetDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            with open(testSetDir + '/appTestingOutput.txt', "w+") as f:
                f.write(out)
                f.write(err)
            
            totalAppEventsFound, totalAppEventsExpected, totalAppGenesFound, totalAppGenesExpected, totalAppEvents = readFiles(testSetDir, 'ApplicationOutput.txt', 'generatorOutput.txt')
            
            if printToConsole:
                print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalAppEventsFound, totalAppEventsExpected, totalAppGenesFound, totalAppGenesExpected, totalAppEvents))
            if totalAppEventsExpected > 0:
                strictAppEventAccuracy = float(totalAppEventsFound)/float(totalAppEventsExpected) * 100.0
            else:
                strictAppEventAccuracy = 0.0
            if totalAppGenesExpected > 0:
                relaxedAppEventAccuracy = float(totalAppGenesFound)/float(totalAppGenesExpected) * 100.0
            else:
                relaxedAppEventAccuracy = 0.0
                
            numEventsAppAveragesList.append(totalAppEvents)
            numEventsGenAveragesList.append(totalAppEventsExpected)
            strictAppAccuracyAveragesList.append(strictAppEventAccuracy)
            relaxedAppAccuracyAveragesList.append(relaxedAppEventAccuracy)
            
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
                if printToConsole:
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
                
                outputEvents(testSetDir + "/orthoAlign.out", testSetDir + "/orthoAlignEvents.out")                
                totalOrthoEventsFound, totalOrthoEventsExpected, totalOrthoGenesFound, totalOrthoGenesExpected, totalOrthoEvents = readFiles(testSetDir, 'orthoAlignEvents.out', 'generatorOutput.txt')
                
                if printToConsole:
                    print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalOrthoEventsFound, totalOrthoEventsExpected, totalOrthoGenesFound, totalOrthoGenesExpected, totalOrthoEvents))
                if totalOrthoEventsExpected > 0:
                    strictOrthoEventAccuracy = float(totalOrthoEventsFound)/float(totalOrthoEventsExpected) * 100.0
                else:
                    strictOrthoEventAccuracy = 0.0
                if totalOrthoGenesExpected > 0:
                    relaxedOrthoEventAccuracy = float(totalOrthoGenesFound)/float(totalOrthoGenesExpected) * 100.0
                else:
                    relaxedOrthoEventAccuracy = 0.0

                strictOrthoAccuracyAveragesList.append(strictOrthoEventAccuracy)
                relaxedOrthoAccuracyAveragesList.append(relaxedOrthoEventAccuracy)
                
                outputEvents(testSetDir + "/duploss.out", testSetDir + "/duplossEvents.out")
                totalDupEventsFound, totalDupEventsExpected, totalDupGenesFound, totalDupGenesExpected, totalDupEvents = readFiles(testSetDir, 'duplossEvents.out', 'generatorOutput.txt')
                
                if printToConsole:
                    print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalDupEventsFound, totalDupEventsExpected, totalDupGenesFound, totalDupGenesExpected, totalDupEvents))
                if totalDupEventsExpected > 0:
                    strictDupEventAccuracy = float(totalDupEventsFound)/float(totalDupEventsExpected) * 100.0
                else:
                    strictDupEventAccuracy = 0.0
                if totalDupGenesExpected > 0:
                    relaxedDupEventAccuracy = float(totalDupGenesFound)/float(totalDupGenesExpected) * 100.0
                else:
                    relaxedDupEventAccuracy = 0.0

                strictDupAccuracyAveragesList.append(strictDupEventAccuracy)
                relaxedDupAccuracyAveragesList.append(relaxedDupEventAccuracy)
                
                if neighbour:
                    p = subprocess.Popen(['python', 'main.py', 'tree2LeafNeighbour.dnd', testSetDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    with open(testSetDir + '/appNeighbourTestingOutput.txt', "w+") as f:
                        f.write(out)
                        f.write(err)
                        
                    totalAppNeighbourEventsFound, totalAppNeighbourEventsExpected, totalAppNeighbourGenesFound, totalAppNeighbourGenesExpected, totalAppNeighbourEvents = readFiles(testSetDir, 'ApplicationNeighbourOutput.txt', 'generatorOutput.txt')
                    
                    if printToConsole:
                        print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalAppNeighbourEventsFound, totalAppNeighbourEventsExpected, totalAppNeighbourGenesFound, totalAppNeighbourGenesExpected, totalAppNeighbourEvents))
                    if totalAppNeighbourEventsExpected > 0:
                        strictAppNeighbourEventAccuracy = float(totalAppNeighbourEventsFound)/float(totalAppNeighbourEventsExpected) * 100.0
                    else:
                        strictAppNeighbourEventAccuracy = 0.0
                    if totalAppNeighbourGenesExpected > 0:
                        relaxedAppNeighbourEventAccuracy = float(totalAppNeighbourGenesFound)/float(totalAppNeighbourGenesExpected) * 100.0
                    else:
                        relaxedAppNeighbourEventAccuracy = 0.0
                        
                    numEventsAppNeighbourAveragesList.append(totalAppNeighbourEvents)
                    strictAppNeighbourAccuracyAveragesList.append(strictAppNeighbourEventAccuracy)
                    relaxedAppNeighbourAccuracyAveragesList.append(relaxedAppNeighbourEventAccuracy)
                    
                    appNeighbourRootFile = testSetDir + "/appNeighbourRoot.txt"
                    with open(appNeighbourRootFile, "r") as f:
                        appNeighbourAncestor = f.readline()
                        
                    appNeighbourRecall, appNeighbourPrecision, appNeighbourfMeasure = compareAnc(appNeighbourAncestor, genAncestor, testSetDir + "/appNeighbour-")
                    appNeighbourFMeasureList.append(appNeighbourfMeasure)
                    
                
        totalEventsAppAveragesList.append(numEventsAppAveragesList)
        totalEventsGenAveragesList.append(numEventsGenAveragesList)
        totalEventsOrthoAveragesList.append(numEventsOrthoAveragesList)
        totalEventsDupAveragesList.append(numEventsDupAveragesList)
        
        totalEventsAppNeighbourAveragesList.append(numEventsAppNeighbourAveragesList)
        
        totalStrictAppAccuracyAveragesList.append(strictAppAccuracyAveragesList)
        totalRelaxedAppAccuracyAveragesList.append(relaxedAppAccuracyAveragesList)
        totalStrictOrthoAccuracyAveragesList.append(strictOrthoAccuracyAveragesList)
        totalRelaxedOrthoAccuracyAveragesList.append(relaxedOrthoAccuracyAveragesList)
        totalStrictDupAccuracyAveragesList.append(strictDupAccuracyAveragesList)
        totalRelaxedDupAccuracyAveragesList.append(relaxedDupAccuracyAveragesList)
        
        totalStrictAppNeighbourAccuracyAveragesList.append(strictAppNeighbourAccuracyAveragesList)
        totalRelaxedAppNeighbourAccuracyAveragesList.append(relaxedAppNeighbourAccuracyAveragesList)
        
        totalAppFMeasureList.append(appFMeasureList)
        totalOrthoFMeasureList.append(orthoFMeasureList)
        totalDupFMeasureList.append(dupFMeasureList)
        
        totalAppNeighbourFMeasureList.append(appNeighbourFMeasureList)
        
        outputData(totalEventsAppAveragesList, testFolder + "appEventsData.txt")
        outputData(totalEventsGenAveragesList, testFolder + "genEventsData.txt")
        outputData(totalEventsOrthoAveragesList, testFolder + "orthoEventsData.txt")
        outputData(totalEventsDupAveragesList, testFolder + "dupEventsData.txt")
        
        outputData(totalEventsAppNeighbourAveragesList, testFolder + "appNeighbourEventsData.txt")
        
        outputData(totalStrictAppAccuracyAveragesList, testFolder + "strictAppAccuracyData.txt")
        outputData(totalRelaxedAppAccuracyAveragesList, testFolder + "relaxedAppAccuracyData.txt")
        outputData(totalStrictOrthoAccuracyAveragesList, testFolder + "strictOrthoAccuracyData.txt")
        outputData(totalRelaxedOrthoAccuracyAveragesList, testFolder + "relaxedOrthoAccuracyData.txt")
        outputData(totalStrictDupAccuracyAveragesList, testFolder + "strictDupAccuracyData.txt")
        outputData(totalRelaxedDupAccuracyAveragesList, testFolder + "relaxedDupAccuracyData.txt")
        
        outputData(totalStrictAppNeighbourAccuracyAveragesList, testFolder + "strictAppNeighbourAccuracyData.txt")
        outputData(totalRelaxedAppNeighbourAccuracyAveragesList, testFolder + "relaxedAppNeighbourAccuracyData.txt")
        
        outputData(totalAppFMeasureList, testFolder + "appFMeasureData.txt")
        outputData(totalOrthoFMeasureList, testFolder + "orthoFMeasureData.txt")
        outputData(totalDupFMeasureList, testFolder + "dupFMeasureData.txt")
        
        outputData(totalAppNeighbourFMeasureList, testFolder + "appNeighbourFMeasureData.txt")
        
        if cherryTree:
            graphData("sAccuracy", totalStrictAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoAccuracyAveragesList, totalAverages4 = totalStrictDupAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourAccuracyAveragesList)
            graphData("rAccuracy", totalRelaxedAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoAccuracyAveragesList, totalAverages4 = totalRelaxedDupAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourAccuracyAveragesList)
            graphData("fMeasure", totalAppFMeasureList, xAxisTitle, xAxis, totalAverages3 = totalOrthoFMeasureList, totalAverages4 = totalDupFMeasureList, totalAverages5 = totalAppNeighbourFMeasureList)
            graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList, totalEventsOrthoAveragesList, totalEventsDupAveragesList, totalEventsAppNeighbourAveragesList)
        else:
            graphData("sAccuracy", totalStrictAppAccuracyAveragesList, xAxisTitle, xAxis)
            graphData("rAccuracy", totalRelaxedAppAccuracyAveragesList, xAxisTitle, xAxis)
            graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList)

        
    if testFolder:
        copy(testFile, testFolder)

def graphData(graphType, totalAverages, xAxisTitle, xAxis, totalAverages2 = None, totalAverages3 = None, totalAverages4 = None, totalAverages5 = None, totalAverages6 = None):
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
    if printToConsole:
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
        if printToConsole:
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
        if printToConsole:
            print totalAverages3
        for averagesList in totalAverages3:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages3.append(average)
        line3, = plt.plot(xAxis, averages3, 'g+-', label='OrthoAlign')
        labels.append(line3)
        
        
    if totalAverages4 is not None:
        averages4 = []
        if printToConsole:
            print totalAverages4
        for averagesList in totalAverages4:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages4.append(average)
        line4, = plt.plot(xAxis, averages4, 'mx-', label='DupLoss')
        labels.append(line4)
        
    if totalAverages5 is not None:
        averages5 = []
        if printToConsole:
            print totalAverages5
        for averagesList in totalAverages5:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages5.append(average)
        line5, = plt.plot(xAxis, averages5, 'o--', label='Application with Neighbour')
        labels.append(line5)
        
    if totalAverages6 is not None:
        averages6 = []
        if printToConsole:
            print totalAverages6
        for averagesList in totalAverages6:
            currentSum = 0.0    
            for average in averagesList:
                currentSum += average
            
            average = currentSum / len(averagesList)
            averages6.append(average)
        line6, = plt.plot(xAxis, averages6, 'mx-', label='DupLoss')
        labels.append(line6)
        
    plt.legend(handles=labels)
    
    if printToConsole:
        print averages
#    plt.show()
    f.savefig(testFolder + graphType + ".pdf", bbox_inches='tight')
    
def outputData(totalAverages, fileName):
    with open(fileName, 'w+') as f:
        for averagesList in totalAverages:
            for average in averagesList:
                f.write("%s " % (average))
            f.write("\n")

main()