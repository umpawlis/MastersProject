from genomeGenerator import generateTests
from CompareResultsService import readFiles
from CompareAncestors import compareAnc
from CompareAncestors import getFirstLineFromFile
from CompareAncestors import cleanUpGenomes
from ListEvents import outputEvents
from shutil import copy
import matplotlib.pyplot as plt
import os
import sys
import subprocess
import datetime
import time

### CONSTANTS ###
ORTHOALIGN_PATH =  "OrthoAlign/OrthoAlign/"; ##I recommend using an absolute path
ORTHOALIGN_EXEC = "Aligning"
DUPLOSS_PATH = "2-SPP/"  ##I recommend using an absolute path
DUPLOSS_EXEC = "duploss"
printToConsole = False

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
    totalEventsOrthoNeighbourAveragesList = []
    
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
    totalStrictOrthoNeighbourAccuracyAveragesList = []
    totalRelaxedOrthoNeighbourAccuracyAveragesList = []
    
    #App total accuracies for different event types
    totalStrictAppDupAccuracyAveragesList = []
    totalRelaxedAppDupAccuracyAveragesList = []
    totalStrictAppLossAccuracyAveragesList = []
    totalRelaxedAppLossAccuracyAveragesList = []
    totalStrictAppInvAccuracyAveragesList = []
    totalRelaxedAppInvAccuracyAveragesList = []
    totalStrictAppTransAccuracyAveragesList = []
    totalRelaxedAppTransAccuracyAveragesList = []
#    totalStrictAppInvTransAccuracyAveragesList = []
#    totalRelaxedAppInvTransAccuracyAveragesList = []
    
    #orthoAlign total accuracies for different event types
    totalStrictOrthoDupAccuracyAveragesList = []
    totalRelaxedOrthoDupAccuracyAveragesList = []
    totalStrictOrthoLossAccuracyAveragesList = []
    totalRelaxedOrthoLossAccuracyAveragesList = []
    totalStrictOrthoInvAccuracyAveragesList = []
    totalRelaxedOrthoInvAccuracyAveragesList = []
    totalStrictOrthoTransAccuracyAveragesList = []
    totalRelaxedOrthoTransAccuracyAveragesList = []
#    totalStrictOrthoInvTransAccuracyAveragesList = []
#    totalRelaxedOrthoInvTransAccuracyAveragesList = []
    
    #Duploss total accuracies for different event types
    totalStrictDupDupAccuracyAveragesList = []
    totalRelaxedDupDupAccuracyAveragesList = []
    totalStrictDupLossAccuracyAveragesList = []
    totalRelaxedDupLossAccuracyAveragesList = []
    totalStrictDupInvAccuracyAveragesList = []
    totalRelaxedDupInvAccuracyAveragesList = []
    totalStrictDupTransAccuracyAveragesList = []
    totalRelaxedDupTransAccuracyAveragesList = []
#    totalStrictDupInvTransAccuracyAveragesList = []
#    totalRelaxedDupInvTransAccuracyAveragesList = []
    
    #AppNeighbour total accuracies for different event types
    totalStrictAppNeighbourDupAccuracyAveragesList = []
    totalRelaxedAppNeighbourDupAccuracyAveragesList = []
    totalStrictAppNeighbourLossAccuracyAveragesList = []
    totalRelaxedAppNeighbourLossAccuracyAveragesList = []
    totalStrictAppNeighbourInvAccuracyAveragesList = []
    totalRelaxedAppNeighbourInvAccuracyAveragesList = []
    totalStrictAppNeighbourTransAccuracyAveragesList = []
    totalRelaxedAppNeighbourTransAccuracyAveragesList = []
#    totalStrictAppNeighbourInvTransAccuracyAveragesList = []
#    totalRelaxedAppNeighbourInvTransAccuracyAveragesList = []
    
    #OrthoNeighbour total accuracies for different event types
    totalStrictOrthoNeighbourDupAccuracyAveragesList = []
    totalRelaxedOrthoNeighbourDupAccuracyAveragesList = []
    totalStrictOrthoNeighbourLossAccuracyAveragesList = []
    totalRelaxedOrthoNeighbourLossAccuracyAveragesList = []
    totalStrictOrthoNeighbourInvAccuracyAveragesList = []
    totalRelaxedOrthoNeighbourInvAccuracyAveragesList = []
    totalStrictOrthoNeighbourTransAccuracyAveragesList = []
    totalRelaxedOrthoNeighbourTransAccuracyAveragesList = []
#    totalStrictOrthoNeighbourInvTransAccuracyAveragesList = []
#    totalRelaxedOrthoNeighbourInvTransAccuracyAveragesList = []
    
    totalAppFMeasureList = []
    totalOrthoFMeasureList = []
    totalDupFMeasureList = []
    
    totalAppNeighbourFMeasureList = []
    totalOrthoNeighbourFMeasureList = []
    
    averageRunTimePerTest = []
    
    for test in tests:
        numEventsAppAveragesList = []
        numEventsGenAveragesList = []
        numEventsOrthoAveragesList = []
        numEventsDupAveragesList = []
        numEventsAppNeighbourAveragesList = []
        numEventsOrthoNeighbourAveragesList = []
        
        strictAppAccuracyAveragesList = []
        relaxedAppAccuracyAveragesList = []
        strictOrthoAccuracyAveragesList = []
        relaxedOrthoAccuracyAveragesList = []
        strictDupAccuracyAveragesList = []
        relaxedDupAccuracyAveragesList = []
        
        strictAppNeighbourAccuracyAveragesList = []
        relaxedAppNeighbourAccuracyAveragesList = []
        strictOrthoNeighbourAccuracyAveragesList = []
        relaxedOrthoNeighbourAccuracyAveragesList = []
        
        #App accuracies for different event types
        strictAppDupAccuracyAveragesList = []
        relaxedAppDupAccuracyAveragesList = []
        strictAppLossAccuracyAveragesList = []
        relaxedAppLossAccuracyAveragesList = []
        strictAppInvAccuracyAveragesList = []
        relaxedAppInvAccuracyAveragesList = []
        strictAppTransAccuracyAveragesList = []
        relaxedAppTransAccuracyAveragesList = []
#        strictAppInvTransAccuracyAveragesList = []
#        relaxedAppInvTransAccuracyAveragesList = []
        
        #orthoAlign accuracies for different event types
        strictOrthoDupAccuracyAveragesList = []
        relaxedOrthoDupAccuracyAveragesList = []
        strictOrthoLossAccuracyAveragesList = []
        relaxedOrthoLossAccuracyAveragesList = []
        strictOrthoInvAccuracyAveragesList = []
        relaxedOrthoInvAccuracyAveragesList = []
        strictOrthoTransAccuracyAveragesList = []
        relaxedOrthoTransAccuracyAveragesList = []
#        strictOrthoInvTransAccuracyAveragesList = []
#        relaxedOrthoInvTransAccuracyAveragesList = []
        
        #duploss accuracies for different event types
        strictDupDupAccuracyAveragesList = []
        relaxedDupDupAccuracyAveragesList = []
        strictDupLossAccuracyAveragesList = []
        relaxedDupLossAccuracyAveragesList = []
        strictDupInvAccuracyAveragesList = []
        relaxedDupInvAccuracyAveragesList = []
        strictDupTransAccuracyAveragesList = []
        relaxedDupTransAccuracyAveragesList = []
#        strictDupInvTransAccuracyAveragesList = []
#        relaxedDupInvTransAccuracyAveragesList = []
        
        #AppNeighbour accuracies for different event types
        strictAppNeighbourDupAccuracyAveragesList = []
        relaxedAppNeighbourDupAccuracyAveragesList = []
        strictAppNeighbourLossAccuracyAveragesList = []
        relaxedAppNeighbourLossAccuracyAveragesList = []
        strictAppNeighbourInvAccuracyAveragesList = []
        relaxedAppNeighbourInvAccuracyAveragesList = []
        strictAppNeighbourTransAccuracyAveragesList = []
        relaxedAppNeighbourTransAccuracyAveragesList = []
#        strictAppNeighbourInvTransAccuracyAveragesList = []
#        relaxedAppNeighbourInvTransAccuracyAveragesList = []
        
        #orthoNeighbour accuracies for different event types
        strictOrthoNeighbourDupAccuracyAveragesList = []
        relaxedOrthoNeighbourDupAccuracyAveragesList = []
        strictOrthoNeighbourLossAccuracyAveragesList = []
        relaxedOrthoNeighbourLossAccuracyAveragesList = []
        strictOrthoNeighbourInvAccuracyAveragesList = []
        relaxedOrthoNeighbourInvAccuracyAveragesList = []
        strictOrthoNeighbourTransAccuracyAveragesList = []
        relaxedOrthoNeighbourTransAccuracyAveragesList = []
#        strictOrthoNeighbourInvTransAccuracyAveragesList = []
#        relaxedOrthoNeighbourInvTransAccuracyAveragesList = []
        
        appFMeasureList = []
        orthoFMeasureList = []
        dupFMeasureList = []
        
        appNeighbourFMeasureList = []
        orthoNeighbourFMeasureList = []
        
        testRunTimes = []
        
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
            startTime = time.time()
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
            
            totalAppEventsFound, totalAppEventsExpected, totalAppGenesFound, totalAppGenesExpected, totalAppEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals = readFiles(testSetDir, 'ApplicationOutput.txt', 'generatorOutput.txt', 'app-')
            strictAppDupEventAccuracy, relaxedAppDupEventAccuracy = calculateAccuracy(duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3])
            strictAppLossEventAccuracy, relaxedAppLossEventAccuracy = calculateAccuracy(lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3])
            strictAppInvEventAccuracy, relaxedAppInvEventAccuracy = calculateAccuracy(inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3])
            strictAppTransEventAccuracy, relaxedAppTransEventAccuracy = calculateAccuracy(transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3])
#            strictAppInvTransEventAccuracy, relaxedAppInvTransEventAccuracy = calculateAccuracy(invertedTranspositionTotals[0], invertedTranspositionTotals[1], invertedTranspositionTotals[2], invertedTranspositionTotals[3])
            
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
            
            strictAppDupAccuracyAveragesList.append(strictAppDupEventAccuracy)
            relaxedAppDupAccuracyAveragesList.append(relaxedAppDupEventAccuracy)
            strictAppLossAccuracyAveragesList.append(strictAppLossEventAccuracy)
            relaxedAppLossAccuracyAveragesList.append(relaxedAppLossEventAccuracy)
            strictAppInvAccuracyAveragesList.append(strictAppInvEventAccuracy)
            relaxedAppInvAccuracyAveragesList.append(relaxedAppInvEventAccuracy)
            strictAppTransAccuracyAveragesList.append(strictAppTransEventAccuracy)
            relaxedAppTransAccuracyAveragesList.append(relaxedAppTransEventAccuracy)
#            strictAppInvTransAccuracyAveragesList.append(strictAppInvTransEventAccuracy)
#            relaxedAppInvTransAccuracyAveragesList.append(relaxedAppInvTransEventAccuracy)
            
            if cherryTree:
                appCommand = baseCommand + testSetDir + '/NC_000001/sequence.txt ' + testSetDir + '/NC_000002/sequence.txt ' + testSetDir
                os.system(appCommand)
                
                duplossOutFile = testSetDir + "/duploss.out"
                orthoAlignOutFile = testSetDir + "/orthoAlign.out"
                appRootFile = testSetDir + "/appRoot.txt"
                if neighbour:
                    genRootFile = testSetDir + "/genAncestor1.txt"
                else:
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
                totalOrthoEventsFound, totalOrthoEventsExpected, totalOrthoGenesFound, totalOrthoGenesExpected, totalOrthoEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals = readFiles(testSetDir, 'orthoAlignEvents.out', 'generatorOutput.txt', 'ortho-')
                strictOrthoDupEventAccuracy, relaxedOrthoDupEventAccuracy = calculateAccuracy(duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3])
                strictOrthoLossEventAccuracy, relaxedOrthoLossEventAccuracy = calculateAccuracy(lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3])
                strictOrthoInvEventAccuracy, relaxedOrthoInvEventAccuracy = calculateAccuracy(inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3])
                strictOrthoTransEventAccuracy, relaxedOrthoTransEventAccuracy = calculateAccuracy(transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3])
                
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
                
                strictOrthoDupAccuracyAveragesList.append(strictOrthoDupEventAccuracy)
                relaxedOrthoDupAccuracyAveragesList.append(relaxedOrthoDupEventAccuracy)
                strictOrthoLossAccuracyAveragesList.append(strictOrthoLossEventAccuracy)
                relaxedOrthoLossAccuracyAveragesList.append(relaxedOrthoLossEventAccuracy)
                strictOrthoInvAccuracyAveragesList.append(strictOrthoInvEventAccuracy)
                relaxedOrthoInvAccuracyAveragesList.append(relaxedOrthoInvEventAccuracy)
                strictOrthoTransAccuracyAveragesList.append(strictOrthoTransEventAccuracy)
                relaxedOrthoTransAccuracyAveragesList.append(relaxedOrthoTransEventAccuracy)
                
                outputEvents(testSetDir + "/duploss.out", testSetDir + "/duplossEvents.out")
                totalDupEventsFound, totalDupEventsExpected, totalDupGenesFound, totalDupGenesExpected, totalDupEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals = readFiles(testSetDir, 'duplossEvents.out', 'generatorOutput.txt', 'dup-')
                strictDupDupEventAccuracy, relaxedDupDupEventAccuracy = calculateAccuracy(duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3])
                strictDupLossEventAccuracy, relaxedDupLossEventAccuracy = calculateAccuracy(lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3])
                strictDupInvEventAccuracy, relaxedDupInvEventAccuracy = calculateAccuracy(inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3])
                strictDupTransEventAccuracy, relaxedDupTransEventAccuracy = calculateAccuracy(transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3])
                
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
                
                strictDupDupAccuracyAveragesList.append(strictDupDupEventAccuracy)
                relaxedDupDupAccuracyAveragesList.append(relaxedDupDupEventAccuracy)
                strictDupLossAccuracyAveragesList.append(strictDupLossEventAccuracy)
                relaxedDupLossAccuracyAveragesList.append(relaxedDupLossEventAccuracy)
                strictDupInvAccuracyAveragesList.append(strictDupInvEventAccuracy)
                relaxedDupInvAccuracyAveragesList.append(relaxedDupInvEventAccuracy)
                strictDupTransAccuracyAveragesList.append(strictDupTransEventAccuracy)
                relaxedDupTransAccuracyAveragesList.append(relaxedDupTransEventAccuracy)
                
                if neighbour:
                    p = subprocess.Popen(['python', 'main.py', 'tree2LeafNeighbour.dnd', testSetDir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    with open(testSetDir + '/appNeighbourTestingOutput.txt', "w+") as f:
                        f.write(out)
                        f.write(err)
                        
                    totalAppNeighbourEventsFound, totalAppNeighbourEventsExpected, totalAppNeighbourGenesFound, totalAppNeighbourGenesExpected, totalAppNeighbourEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals = readFiles(testSetDir, 'ApplicationNeighbourOutput.txt', 'generatorOutput.txt', 'appNeighbour-')
                    strictAppNeighbourDupEventAccuracy, relaxedAppNeighbourDupEventAccuracy = calculateAccuracy(duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3])
                    strictAppNeighbourLossEventAccuracy, relaxedAppNeighbourLossEventAccuracy = calculateAccuracy(lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3])
                    strictAppNeighbourInvEventAccuracy, relaxedAppNeighbourInvEventAccuracy = calculateAccuracy(inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3])
                    strictAppNeighbourTransEventAccuracy, relaxedAppNeighbourTransEventAccuracy = calculateAccuracy(transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3])
#                    strictAppNeighbourInvTransEventAccuracy, relaxedAppNeighbourInvTransEventAccuracy = calculateAccuracy(invertedTranspositionTotals[0], invertedTranspositionTotals[1], invertedTranspositionTotals[2], invertedTranspositionTotals[3])
                    
                    
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
                    
                    strictAppNeighbourDupAccuracyAveragesList.append(strictAppNeighbourDupEventAccuracy)
                    relaxedAppNeighbourDupAccuracyAveragesList.append(relaxedAppNeighbourDupEventAccuracy)
                    strictAppNeighbourLossAccuracyAveragesList.append(strictAppNeighbourLossEventAccuracy)
                    relaxedAppNeighbourLossAccuracyAveragesList.append(relaxedAppNeighbourLossEventAccuracy)
                    strictAppNeighbourInvAccuracyAveragesList.append(strictAppNeighbourInvEventAccuracy)
                    relaxedAppNeighbourInvAccuracyAveragesList.append(relaxedAppNeighbourInvEventAccuracy)
                    strictAppNeighbourTransAccuracyAveragesList.append(strictAppNeighbourTransEventAccuracy)
                    relaxedAppNeighbourTransAccuracyAveragesList.append(relaxedAppNeighbourTransEventAccuracy)
#                    strictAppNeighbourInvTransAccuracyAveragesList.append(strictAppNeighbourInvTransEventAccuracy)
#                    relaxedAppNeighbourInvTransAccuracyAveragesList.append(relaxedAppNeighbourInvTransEventAccuracy)
                    
                    appNeighbourRootFile = testSetDir + "/appNeighbourRoot.txt"
                    with open(appNeighbourRootFile, "r") as f:
                        appNeighbourAncestor = f.readline()
                        
                    appNeighbourRecall, appNeighbourPrecision, appNeighbourfMeasure = compareAnc(appNeighbourAncestor, genAncestor, testSetDir + "/appNeighbour-")
                    appNeighbourFMeasureList.append(appNeighbourfMeasure)
                    
                    genome1 = getFirstLineFromFile(testSetDir + '/NC_000001/sequence.txt ')
                    genome2 = getFirstLineFromFile(testSetDir + '/NC_000002/sequence.txt ')
                    genome3 = getFirstLineFromFile(testSetDir + '/NC_000003/sequence.txt ')
                    
                    genome1 = cleanUpGenomes(genome1)
                    genome2 = cleanUpGenomes(genome2)
                    genome3 = cleanUpGenomes(genome3)
                    
                    #Running Duploss with neighbour
                    command = "java -classpath " + ORTHOALIGN_PATH + " " + ORTHOALIGN_EXEC + " -dt " + genome1 + " " + genome2 + " " + genome3 + " > " + testSetDir + "/orthoAlignNeighbour.out"
                    print command
                    os.system(command)
                    
                    orthoAlignNeighbourOutFile = testSetDir + "/orthoAlignNeighbour.out"
                    
                    with open(orthoAlignNeighbourOutFile, "r") as f:
                        line = f.readline()
                        while line:
                            splitted = line.split("ost = ")
                            if len(splitted) > 1:
                                orthoNeighbourCost = float(splitted[1])
                                line = f.readline()
                            if line.strip() == ">Ancestor:":
                                orthoNeighbourAncestor = f.readline().strip()
                                break
                            line = f.readline()
                            
                    numEventsOrthoNeighbourAveragesList.append(orthoNeighbourCost)
                    orthoNeighbourRecall, orthoNeighbourPrecision, orthoNeighbourfMeasure = compareAnc(orthoNeighbourAncestor, genAncestor, testSetDir + "/orthoNeighbour-")
                    orthoNeighbourFMeasureList.append(orthoNeighbourfMeasure)
                    
                    outputEvents(testSetDir + "/orthoAlignNeighbour.out", testSetDir + "/orthoAlignNeighbourEvents.out")                
                    totalOrthoNeighbourEventsFound, totalOrthoNeighbourEventsExpected, totalOrthoNeighbourGenesFound, totalOrthoNeighbourGenesExpected, totalOrthoNeighbourEvents, duplicationTotals, lossTotals, inversionTotals, transpositionTotals = readFiles(testSetDir, 'orthoAlignNeighbourEvents.out', 'generatorOutput.txt', 'orthoNeighbour-')
                    strictOrthoNeighbourDupEventAccuracy, relaxedOrthoNeighbourDupEventAccuracy = calculateAccuracy(duplicationTotals[0], duplicationTotals[1], duplicationTotals[2], duplicationTotals[3])
                    strictOrthoNeighbourLossEventAccuracy, relaxedOrthoNeighbourLossEventAccuracy = calculateAccuracy(lossTotals[0], lossTotals[1], lossTotals[2], lossTotals[3])
                    strictOrthoNeighbourInvEventAccuracy, relaxedOrthoNeighbourInvEventAccuracy = calculateAccuracy(inversionTotals[0], inversionTotals[1], inversionTotals[2], inversionTotals[3])
                    strictOrthoNeighbourTransEventAccuracy, relaxedOrthoNeighbourTransEventAccuracy = calculateAccuracy(transpositionTotals[0], transpositionTotals[1], transpositionTotals[2], transpositionTotals[3])
                    
                    if printToConsole:
                        print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s Total App Events: %s' % (totalOrthoNeighbourEventsFound, totalOrthoNeighbourEventsExpected, totalOrthoNeighbourGenesFound, totalOrthoNeighbourGenesExpected, totalOrthoNeighbourEvents))
                    if totalOrthoNeighbourEventsExpected > 0:
                        strictOrthoNeighbourEventAccuracy = float(totalOrthoNeighbourEventsFound)/float(totalOrthoNeighbourEventsExpected) * 100.0
                    else:
                        strictOrthoNeighbourEventAccuracy = 0.0
                    if totalOrthoNeighbourGenesExpected > 0:
                        relaxedOrthoNeighbourEventAccuracy = float(totalOrthoNeighbourGenesFound)/float(totalOrthoNeighbourGenesExpected) * 100.0
                    else:
                        relaxedOrthoNeighbourEventAccuracy = 0.0
    
                    strictOrthoNeighbourAccuracyAveragesList.append(strictOrthoNeighbourEventAccuracy)
                    relaxedOrthoNeighbourAccuracyAveragesList.append(relaxedOrthoNeighbourEventAccuracy)
                    
                    strictOrthoNeighbourDupAccuracyAveragesList.append(strictOrthoNeighbourDupEventAccuracy)
                    relaxedOrthoNeighbourDupAccuracyAveragesList.append(relaxedOrthoNeighbourDupEventAccuracy)
                    strictOrthoNeighbourLossAccuracyAveragesList.append(strictOrthoNeighbourLossEventAccuracy)
                    relaxedOrthoNeighbourLossAccuracyAveragesList.append(relaxedOrthoNeighbourLossEventAccuracy)
                    strictOrthoNeighbourInvAccuracyAveragesList.append(strictOrthoNeighbourInvEventAccuracy)
                    relaxedOrthoNeighbourInvAccuracyAveragesList.append(relaxedOrthoNeighbourInvEventAccuracy)
                    strictOrthoNeighbourTransAccuracyAveragesList.append(strictOrthoNeighbourTransEventAccuracy)
                    relaxedOrthoNeighbourTransAccuracyAveragesList.append(relaxedOrthoNeighbourTransEventAccuracy)
            runTime = time.time() - startTime
            testRunTimes.append(runTime)
            
        averageRunTimePerTest.append(testRunTimes)
        printAverages(averageRunTimePerTest)
        
        totalEventsAppAveragesList.append(numEventsAppAveragesList)
        totalEventsGenAveragesList.append(numEventsGenAveragesList)
        totalEventsOrthoAveragesList.append(numEventsOrthoAveragesList)
        totalEventsDupAveragesList.append(numEventsDupAveragesList)
        
        totalEventsAppNeighbourAveragesList.append(numEventsAppNeighbourAveragesList)
        totalEventsOrthoNeighbourAveragesList.append(numEventsOrthoNeighbourAveragesList)
        
        totalStrictAppAccuracyAveragesList.append(strictAppAccuracyAveragesList)
        totalRelaxedAppAccuracyAveragesList.append(relaxedAppAccuracyAveragesList)
        totalStrictOrthoAccuracyAveragesList.append(strictOrthoAccuracyAveragesList)
        totalRelaxedOrthoAccuracyAveragesList.append(relaxedOrthoAccuracyAveragesList)
        totalStrictDupAccuracyAveragesList.append(strictDupAccuracyAveragesList)
        totalRelaxedDupAccuracyAveragesList.append(relaxedDupAccuracyAveragesList)
        
        totalStrictAppNeighbourAccuracyAveragesList.append(strictAppNeighbourAccuracyAveragesList)
        totalRelaxedAppNeighbourAccuracyAveragesList.append(relaxedAppNeighbourAccuracyAveragesList)
        totalStrictOrthoNeighbourAccuracyAveragesList.append(strictOrthoNeighbourAccuracyAveragesList)
        totalRelaxedOrthoNeighbourAccuracyAveragesList.append(relaxedOrthoNeighbourAccuracyAveragesList)
        
        totalStrictAppDupAccuracyAveragesList.append(strictAppDupAccuracyAveragesList)
        totalRelaxedAppDupAccuracyAveragesList.append(relaxedAppDupAccuracyAveragesList)
        totalStrictAppLossAccuracyAveragesList.append(strictAppLossAccuracyAveragesList)
        totalRelaxedAppLossAccuracyAveragesList.append(relaxedAppLossAccuracyAveragesList)
        totalStrictAppInvAccuracyAveragesList.append(strictAppInvAccuracyAveragesList)
        totalRelaxedAppInvAccuracyAveragesList.append(relaxedAppInvAccuracyAveragesList)
        totalStrictAppTransAccuracyAveragesList.append(strictAppTransAccuracyAveragesList)
        totalRelaxedAppTransAccuracyAveragesList.append(relaxedAppTransAccuracyAveragesList)
#        totalStrictAppInvTransAccuracyAveragesList.append(strictAppInvTransAccuracyAveragesList)
#        totalRelaxedAppInvTransAccuracyAveragesList.append(relaxedAppInvTransAccuracyAveragesList)
        
        totalStrictOrthoDupAccuracyAveragesList.append(strictOrthoDupAccuracyAveragesList)
        totalRelaxedOrthoDupAccuracyAveragesList.append(relaxedOrthoDupAccuracyAveragesList)
        totalStrictOrthoLossAccuracyAveragesList.append(strictOrthoLossAccuracyAveragesList)
        totalRelaxedOrthoLossAccuracyAveragesList.append(relaxedOrthoLossAccuracyAveragesList)
        totalStrictOrthoInvAccuracyAveragesList.append(strictOrthoInvAccuracyAveragesList)
        totalRelaxedOrthoInvAccuracyAveragesList.append(relaxedOrthoInvAccuracyAveragesList)
        totalStrictOrthoTransAccuracyAveragesList.append(strictOrthoTransAccuracyAveragesList)
        totalRelaxedOrthoTransAccuracyAveragesList.append(relaxedOrthoTransAccuracyAveragesList)
        
        totalStrictDupDupAccuracyAveragesList.append(strictDupDupAccuracyAveragesList)
        totalRelaxedDupDupAccuracyAveragesList.append(relaxedDupDupAccuracyAveragesList)
        totalStrictDupLossAccuracyAveragesList.append(strictDupLossAccuracyAveragesList)
        totalRelaxedDupLossAccuracyAveragesList.append(relaxedDupLossAccuracyAveragesList)
        totalStrictDupInvAccuracyAveragesList.append(strictDupInvAccuracyAveragesList)
        totalRelaxedDupInvAccuracyAveragesList.append(relaxedDupInvAccuracyAveragesList)
        totalStrictDupTransAccuracyAveragesList.append(strictDupTransAccuracyAveragesList)
        totalRelaxedDupTransAccuracyAveragesList.append(relaxedDupTransAccuracyAveragesList)
        
        totalStrictAppNeighbourDupAccuracyAveragesList.append(strictAppNeighbourDupAccuracyAveragesList)
        totalRelaxedAppNeighbourDupAccuracyAveragesList.append(relaxedAppNeighbourDupAccuracyAveragesList)
        totalStrictAppNeighbourLossAccuracyAveragesList.append(strictAppNeighbourLossAccuracyAveragesList)
        totalRelaxedAppNeighbourLossAccuracyAveragesList.append(relaxedAppNeighbourLossAccuracyAveragesList)
        totalStrictAppNeighbourInvAccuracyAveragesList.append(strictAppNeighbourInvAccuracyAveragesList)
        totalRelaxedAppNeighbourInvAccuracyAveragesList.append(relaxedAppNeighbourInvAccuracyAveragesList)
        totalStrictAppNeighbourTransAccuracyAveragesList.append(strictAppNeighbourTransAccuracyAveragesList)
        totalRelaxedAppNeighbourTransAccuracyAveragesList.append(relaxedAppNeighbourTransAccuracyAveragesList)
#        totalStrictAppNeighbourInvTransAccuracyAveragesList.append(strictAppNeighbourInvTransAccuracyAveragesList)
#        totalRelaxedAppNeighbourInvTransAccuracyAveragesList.append(relaxedAppNeighbourInvTransAccuracyAveragesList)
        
        totalStrictOrthoNeighbourDupAccuracyAveragesList.append(strictOrthoNeighbourDupAccuracyAveragesList)
        totalRelaxedOrthoNeighbourDupAccuracyAveragesList.append(relaxedOrthoNeighbourDupAccuracyAveragesList)
        totalStrictOrthoNeighbourLossAccuracyAveragesList.append(strictOrthoNeighbourLossAccuracyAveragesList)
        totalRelaxedOrthoNeighbourLossAccuracyAveragesList.append(relaxedOrthoNeighbourLossAccuracyAveragesList)
        totalStrictOrthoNeighbourInvAccuracyAveragesList.append(strictOrthoNeighbourInvAccuracyAveragesList)
        totalRelaxedOrthoNeighbourInvAccuracyAveragesList.append(relaxedOrthoNeighbourInvAccuracyAveragesList)
        totalStrictOrthoNeighbourTransAccuracyAveragesList.append(strictOrthoNeighbourTransAccuracyAveragesList)
        totalRelaxedOrthoNeighbourTransAccuracyAveragesList.append(relaxedOrthoNeighbourTransAccuracyAveragesList)
        
        totalAppFMeasureList.append(appFMeasureList)
        totalOrthoFMeasureList.append(orthoFMeasureList)
        totalDupFMeasureList.append(dupFMeasureList)
        
        totalAppNeighbourFMeasureList.append(appNeighbourFMeasureList)
        totalOrthoNeighbourFMeasureList.append(orthoFMeasureList)
        
        outputData(totalEventsAppAveragesList, testFolder + "appEventsData.txt")
        outputData(totalEventsGenAveragesList, testFolder + "genEventsData.txt")
        outputData(totalEventsOrthoAveragesList, testFolder + "orthoEventsData.txt")
        outputData(totalEventsDupAveragesList, testFolder + "dupEventsData.txt")
        
        outputData(totalEventsAppNeighbourAveragesList, testFolder + "appNeighbourEventsData.txt")
        outputData(totalEventsOrthoNeighbourAveragesList, testFolder + "orthoNeighbourEventsData.txt")
        
        outputData(totalStrictAppAccuracyAveragesList, testFolder + "strictAppAccuracyData.txt")
        outputData(totalRelaxedAppAccuracyAveragesList, testFolder + "relaxedAppAccuracyData.txt")
        outputData(totalStrictOrthoAccuracyAveragesList, testFolder + "strictOrthoAccuracyData.txt")
        outputData(totalRelaxedOrthoAccuracyAveragesList, testFolder + "relaxedOrthoAccuracyData.txt")
        outputData(totalStrictDupAccuracyAveragesList, testFolder + "strictDupAccuracyData.txt")
        outputData(totalRelaxedDupAccuracyAveragesList, testFolder + "relaxedDupAccuracyData.txt")
        
        outputData(totalStrictAppNeighbourAccuracyAveragesList, testFolder + "strictAppNeighbourAccuracyData.txt")
        outputData(totalRelaxedAppNeighbourAccuracyAveragesList, testFolder + "relaxedAppNeighbourAccuracyData.txt")
        outputData(totalStrictOrthoNeighbourAccuracyAveragesList, testFolder + "strictOrthoNeighbourAccuracyData.txt")
        outputData(totalRelaxedOrthoNeighbourAccuracyAveragesList, testFolder + "relaxedOrthoNeighbourAccuracyData.txt")
        
        outputData(totalAppFMeasureList, testFolder + "appFMeasureData.txt")
        outputData(totalOrthoFMeasureList, testFolder + "orthoFMeasureData.txt")
        outputData(totalDupFMeasureList, testFolder + "dupFMeasureData.txt")
        
        outputData(totalAppNeighbourFMeasureList, testFolder + "appNeighbourFMeasureData.txt")
        outputData(totalOrthoNeighbourFMeasureList, testFolder + "orthoNeighbourFMeasureData.txt")
        
        if cherryTree:
            if neighbour: 
                graphData("sAccuracy", totalStrictAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoAccuracyAveragesList, totalAverages4 = totalStrictDupAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourAccuracyAveragesList, totalAverages6 = totalStrictOrthoNeighbourAccuracyAveragesList)
                graphData("rAccuracy", totalRelaxedAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoAccuracyAveragesList, totalAverages4 = totalRelaxedDupAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourAccuracyAveragesList, totalAverages6 = totalRelaxedOrthoNeighbourAccuracyAveragesList)
                graphData("fMeasure", totalAppFMeasureList, xAxisTitle, xAxis, totalAverages3 = totalOrthoFMeasureList, totalAverages4 = totalDupFMeasureList, totalAverages5 = totalAppNeighbourFMeasureList, totalAverages6 = totalOrthoNeighbourFMeasureList)
                graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList, totalEventsOrthoAveragesList, totalEventsDupAveragesList, totalEventsAppNeighbourAveragesList, totalEventsOrthoNeighbourAveragesList)
                
                graphData("sAccuracyDup", totalStrictAppDupAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoDupAccuracyAveragesList, totalAverages4 = totalStrictDupDupAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourDupAccuracyAveragesList, totalAverages6 = totalStrictOrthoNeighbourDupAccuracyAveragesList)
                graphData("rAccuracyDup", totalRelaxedAppDupAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoDupAccuracyAveragesList, totalAverages4 = totalRelaxedDupDupAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourDupAccuracyAveragesList, totalAverages6 = totalRelaxedOrthoNeighbourDupAccuracyAveragesList)
                graphData("sAccuracyLoss", totalStrictAppLossAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoLossAccuracyAveragesList, totalAverages4 = totalStrictDupLossAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourLossAccuracyAveragesList, totalAverages6 = totalStrictOrthoNeighbourLossAccuracyAveragesList)
                graphData("rAccuracyLoss", totalRelaxedAppLossAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoLossAccuracyAveragesList, totalAverages4 = totalRelaxedDupLossAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourLossAccuracyAveragesList, totalAverages6 = totalRelaxedOrthoNeighbourLossAccuracyAveragesList)
                graphData("sAccuracyInv", totalStrictAppInvAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoInvAccuracyAveragesList, totalAverages4 = totalStrictDupInvAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourInvAccuracyAveragesList, totalAverages6 = totalStrictOrthoNeighbourInvAccuracyAveragesList)
                graphData("rAccuracyInv", totalRelaxedAppInvAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoInvAccuracyAveragesList, totalAverages4 = totalRelaxedDupInvAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourInvAccuracyAveragesList, totalAverages6 = totalRelaxedOrthoNeighbourInvAccuracyAveragesList)
                graphData("sAccuracyTrans", totalStrictAppTransAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoTransAccuracyAveragesList, totalAverages4 = totalStrictDupTransAccuracyAveragesList, totalAverages5 = totalStrictAppNeighbourTransAccuracyAveragesList, totalAverages6 = totalStrictOrthoNeighbourTransAccuracyAveragesList)
                graphData("rAccuracyTrans", totalRelaxedAppTransAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoTransAccuracyAveragesList, totalAverages4 = totalRelaxedDupTransAccuracyAveragesList, totalAverages5 = totalRelaxedAppNeighbourTransAccuracyAveragesList, totalAverages6 = totalRelaxedOrthoNeighbourTransAccuracyAveragesList)
            else:
                graphData("sAccuracy", totalStrictAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalStrictOrthoAccuracyAveragesList, totalAverages4 = totalStrictDupAccuracyAveragesList)
                graphData("rAccuracy", totalRelaxedAppAccuracyAveragesList, xAxisTitle, xAxis, totalAverages3 = totalRelaxedOrthoAccuracyAveragesList, totalAverages4 = totalRelaxedDupAccuracyAveragesList)
                graphData("fMeasure", totalAppFMeasureList, xAxisTitle, xAxis, totalAverages3 = totalOrthoFMeasureList, totalAverages4 = totalDupFMeasureList)
                graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList, totalEventsOrthoAveragesList, totalEventsDupAveragesList)
        else:
            graphData("sAccuracy", totalStrictAppAccuracyAveragesList, xAxisTitle, xAxis)
            graphData("rAccuracy", totalRelaxedAppAccuracyAveragesList, xAxisTitle, xAxis)
            graphData("Events", totalEventsAppAveragesList, xAxisTitle, xAxis, totalEventsGenAveragesList)

        
    if testFolder:
        copy(testFile, testFolder)
        
def calculateAccuracy(totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected):
    if printToConsole:
        print('Events Found: %s Events Expected: %s Genes Found: %s Genes Expected: %s' % (totalEventsFound, totalEventsExpected, totalGenesFound, totalGenesExpected))
    if totalEventsExpected > 0:
        strictEventAccuracy = float(totalEventsFound)/float(totalEventsExpected) * 100.0
    else:
        strictEventAccuracy = 0.0
    if totalGenesExpected > 0:
        relaxedEventAccuracy = float(totalGenesFound)/float(totalGenesExpected) * 100.0
    else:
        relaxedEventAccuracy = 0.0
        
    return strictEventAccuracy, relaxedEventAccuracy

def printAverages(AveragesPerTest):
    runTimeSum = 0
    
    with open(testFolder + "runtimeAverages.txt", "w+") as f:
        for testRuntimes in AveragesPerTest:
            for runtime in testRuntimes:
                runTimeSum += runtime
                f.write(runtime + " ")
                
            average = runTimeSum / len(testRuntimes)
            f.write("\n")
            f.write(str(average) + "\n")
        

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
        
    elif graphType == "sAccuracyDup":
        title = "Average Strict Duplication Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracyDup":
        title = "Average Relaxed Duplication Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
    elif graphType == "sAccuracyLoss":
        title = "Average Strict Loss Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracyLoss":
        title = "Average Relaxed Loss Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
    elif graphType == "sAccuracyInv":
        title = "Average Strict Inversion Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracyInv":
        title = "Average Relaxed Inversion Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
    elif graphType == "sAccuracyTrans":
        title = "Average Strict Transposition Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracyTrans":
        title = "Average Relaxed Transposition Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
    elif graphType == "sAccuracyInvTrans":
        title = "Average Strict Inverted Transposition Accuracy"
        yAxisTitle = "Accuracy Percentage"
    elif graphType == "rAccuracyInvTrans":
        title = "Average Relaxed Inverted Transposition Accuracy"
        yAxisTitle = "Accuracy Percentage"
        
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
        line6, = plt.plot(xAxis, averages6, 'yx--', label='OrthoAlign with Neighbour')
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