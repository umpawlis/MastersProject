import re
import sys

### CONSTANTS ###



### METHODS ###
def parseOutputFile(filename):
	"""
	Parses either a duploss or an OrthoAlign output file.
	Returns arrays representing all the genes lost, duplicated, inverted, transposed and substituted.
	Also returns the position of the terminus in X and Y, and the position of the last gene in X.
	These positions can be used to calculate the corresponding positions in our application and generator.
	Know that in duploss and OrthoAlign, [o] and [t] count as genes so they are assigned an index (need to correct for that) -->
	if the gene is after [o] but before [t], do -1. If the gene is after [t], do -2.
	You can get a corresponding position in genome Y (the second one) by subtracting the index of the gene with posLastGeneX, and
	correcting for the indexes used by [o] and [t] (if necessary)
	"""
	try:
		f = open(filename, 'r')
	except IOError as e:
		sys.exit("Error opening file "+filename+': '+str(e))
	
	#initializing some lists of events
	dels = []
	dups = []
	invs = []
	transpos = []
	invTranspos = [] #new!
	subs = []
	
	#Keeping track of some positions --> will be useful to convert positions
	posTerminusX = -1  #Not used anymore
	posLastGeneX = -1
	posTerminusY = -1  #Not used anymore
	
	#Keeping track of the genomes --> mapping index to gene
	indexToGeneDict = {}
	
	startListing = False  #start listing events
	line = f.readline()
	while line:
		if not startListing and "========" in line:
			startListing = True
		
		elif startListing and line.strip():  #OrthoAlign has an empty line before the start of the alignment...
			splitted = line.split(" --- ")
			if len(splitted) < 2: #it means we're done parsing the alignment
				break
			
			#if "[t]" in splitted[0]:  #OrthoAlign does not always match [t]... so we have to treat them separately
			#	posTerminusX = int(splitted[0].split("_")[1])
			#	indexToGeneDict[posTerminusX] = "< t >";
			#if ". " not in splitted[1] and "[t]" in splitted[1]:  #OrthoAlign can have an event associated with [t]...
			#	posTerminusY = int(splitted[1].split("_")[1].split()[0])
			#	indexToGeneDict[posTerminusY] = "< t >";
			
			#else:
			if " ." not in splitted[1]:
				geneSplit = splitted[1].split()[0].split("_")
				geneIndex = int(geneSplit[len(geneSplit)-1])
				geneName = geneSplit[0].replace("-", "").replace("[t]", "< t >")
				indexToGeneDict[geneIndex] = geneName  #rRNA genes don't have an anticodon
				if len(geneSplit) > 2:  #for tRNA genes with anticodons...
					indexToGeneDict[geneIndex] += "_" + geneSplit[1]  #add the anticodon
				
			if "." not in splitted[0]:
				geneSplit = splitted[0].replace(" ", "").split("_")
				geneIndex = int(geneSplit[len(geneSplit)-1])
				geneName = geneSplit[0].replace("-", "").replace("[t]", "< t >")
				indexToGeneDict[geneIndex] = geneName   #rRNA genes don't have an anticodon
				if len(geneSplit) > 2:  #for tRNA genes with anticodons...
					indexToGeneDict[geneIndex] += "_" + geneSplit[1]  #add the anticodon
					
				posLastGeneX = geneIndex  #updating at every step
				
			operationSplit = splitted[1].split()[1:]  #skipping the gene, getting the rest (operation description)
			if len(operationSplit) > 0:
				if operationSplit[0] == "Substitution":
					#subs.append(splitted[0].replace(" ", "") + "<->" + splitted[1].split("Substitution")[0].replace(" ", ""))
					#now splitting them into two events --> not anymore, I realized they always assign it to genome Y
					#subs.append(splitted[0].replace(" ", ""))  #gene in genome X is always considered to be the ancestral one
					subs.append(splitted[1].split("Substitution")[0].replace(" ", ""))
				else:
					genesAffected = operationSplit[2]
					startEndSplit = operationSplit[2].split("...")
					if len(startEndSplit) > 1:
						if(startEndSplit[0] == startEndSplit[1]):  #for OrthoAlign --> means only one gene affected
							genesAffected = startEndSplit[0]
					if operationSplit[0] == "Loss":
						dels.append(genesAffected)
					elif operationSplit[0] == "Duplication" or operationSplit[0] == "Inverted": #Inverted is for Inverted duplication
						if "to" in operationSplit:  #that means it's a duploss file
							if operationSplit[2] == operationSplit[4]: #only one gene affected
								genesAffected = operationSplit[2]
							else:
								genesAffected = operationSplit[2] + "..." + operationSplit[4]  #transforming into OrthoAlign format
						if operationSplit[0] == "Inverted":
							if "to" in operationSplit:  #that means it's a duploss file
								if operationSplit[3] == operationSplit[5]: #only one gene affected
									genesAffected = operationSplit[3]
								else:
									genesAffected = operationSplit[3] + "..." + operationSplit[5]  #transforming into OrthoAlign format
							else: #OrthoAlign file
								genesAffected = operationSplit[3]  #Inverted duplication of [genesAffected]
								startEndSplit = operationSplit[2].split("...")
								if len(startEndSplit) > 1:
									if(startEndSplit[0] == startEndSplit[1]):  #for OrthoAlign --> means only one gene affected
										genesAffected = startEndSplit[0]
						dups.append(genesAffected)
					elif operationSplit[0] == "Inversion":
						#invs.append(genesAffected)  --> not anymore, I realized they always assign it to genome Y
						otherGenesAffected = operationSplit[4]
						startEndSplit = operationSplit[4].split("...")
						if len(startEndSplit) > 1:
							if(startEndSplit[0] == startEndSplit[1]):  #for OrthoAlign --> means only one gene affected
								otherGenesAffected = startEndSplit[0]
						invs.append(otherGenesAffected)
					elif operationSplit[0] == "transposition":
						if "." in splitted[0]: #Now adding the transposition only on genome Y, because I realized that's where they always put it
							if ("-" in genesAffected and "-" in operationSplit[4]) or ("-" not in genesAffected and "-" not in operationSplit[4]):  #regular transpo
								transpos.append(genesAffected)
							else: #inverted transpo
								invTranspos.append(genesAffected)
					else:
						print "BIG PROBLEM: operation " + operationSplit[0] + " not recognized"
			
		line = f.readline()
		
	#the prints are for testing purposes
	#print dels
	#print dups
	#print invs
	#print transpos
	#print subs
	#print posTerminusX
	#print posLastGeneX
	#print posTerminusY
	#print indexToGeneDict
	
	return dels, dups, invs, transpos, invTranspos, subs, posTerminusX, posLastGeneX, posTerminusY, indexToGeneDict
	
	
def outputEvents(algoOutput, outputFile):

	f = open(outputFile, 'w')
	
	dels, dups, invs, transpos, invTranspos, subs, posTerminusX, posLastGeneX, posTerminusY, indexToGeneDict = parseOutputFile(algoOutput)
	
	#Header stuff
	eventsStrX = "#tree\n"
	eventsStrX += "Strain:NC_000001\n"
	eventsStrY = "Strain:NC_000002\n"
	eventsStrX += "Codon Mismatch:\n"  #those would not be identified by OrthoAlign nor duploss
	eventsStrY += "Codon Mismatch:\n"  #those would not be identified by OrthoAlign nor duploss
	
	### Substitutions ###
	eventsStrX += "Substitution:"
	eventsStrY += "Substitution:"
	
	for sub in subs: #including them on both sides (both genomes), since OrthoAlign does not specify where it occurred
		#ind = int(re.search("[0-9]+$", sub).group(0))
		splitted = sub.split('_')
		ind = int(splitted[len(splitted)-1])
		s = indexToGeneDict[ind]#splitted[0]
		#if len(splitted) > 2:
		#	s += '_' + splitted[1]  #anticodon
		s += ' ' + str(convertIndex(ind, posLastGeneX)) + ';'
		if ind <= posLastGeneX:
			eventsStrX += s
		else:
			eventsStrY += s
	eventsStrX += "\n"
	eventsStrY += "\n"

	### Duplications ###
	eventsStrX += "Duplication:"
	eventsStrY += "Duplication:"
	
	for dup in dups:
		s, pos = getEventStr(dup, indexToGeneDict, posLastGeneX)
		if pos > posLastGeneX:
			eventsStrY += s
		else: eventsStrX += s
		
	eventsStrX += "\n"
	eventsStrY += "\n"
	
	### Deletions ###
	eventsStrX += "Deletion:"
	eventsStrY += "Deletion:"
	
	for loss in dels:
		s, pos = getEventStr(loss, indexToGeneDict, posLastGeneX)
		if pos < posLastGeneX:  #--> it's the opposite for deletions (if the pos is in Y, the del is in X)
			eventsStrY += s
		else: eventsStrX += s
		
	eventsStrX += "\n"
	eventsStrY += "\n"
	
	### Inversions ###
	eventsStrX += "Inversion:"
	eventsStrY += "Inversion:"
	
	for inv in invs:
		s, pos = getEventStr(inv, indexToGeneDict, posLastGeneX, True)
		if pos > posLastGeneX:
			eventsStrY += s
		else: eventsStrX += s
		
	eventsStrX += "\n"
	eventsStrY += "\n"
	
	
	### Transpositions ###
	eventsStrX += "Transposition:"
	eventsStrY += "Transposition:"
	
	for transpo in transpos:
		s, pos = getEventStr(transpo, indexToGeneDict, posLastGeneX, True)
		if pos > posLastGeneX:
			eventsStrY += s
		else: eventsStrX += s
		
	eventsStrX += "\n"
	eventsStrY += "\n"
	
	### Inverted transpositions ###
	eventsStrX += "Inverted Transposition:"
	eventsStrY += "Inverted Transposition:"
	
	for invTranspo in invTranspos:
		s, pos = getEventStr(invTranspo, indexToGeneDict, posLastGeneX, True)
		if pos > posLastGeneX:
			eventsStrY += s
		else: eventsStrX += s
		
	eventsStrX += "\n"
	eventsStrY += "\n"
	
	try:
		f.write(eventsStrX + eventsStrY)
		f.close();
	
	except IOError as e:
		sys.exit("Error writing to file "+filename+': '+e)
	
	
def getEventStr(event, indexToGeneDict, posLastGeneX, rearr = False): #rearr if it's an inv or transpo
	splitted = event.split("...")
	if len(splitted) > 1:
		startPos = int(re.search("[0-9]+$", splitted[0]).group(0))
		endPos = int(re.search("[0-9]+$", splitted[1]).group(0))
			
		eventStr = ""
		for pos in range(startPos, endPos):
			eventStr += indexToGeneDict[pos] + " " + str(convertIndex(pos, posLastGeneX)) + ", "
		eventStr += indexToGeneDict[endPos] + " " + str(convertIndex(endPos, posLastGeneX))
	else:
		splitted = event.split('_')
		ind = int(splitted[len(splitted)-1])
		#eventStr = splitted[0]
		#if len(splitted) > 2:
		#	eventStr += '_' + splitted[1]  #anticodon
		#eventStr += ' ' + str(ind) + ';'
		endPos = ind
		eventStr = indexToGeneDict[ind] + ' ' + str(convertIndex(ind, posLastGeneX))
	if rearr:
		eventStr += ';|'
	else:
		eventStr += ';'
	return eventStr, endPos
	
def convertIndex(index, posLastGeneX):
	if index > posLastGeneX:
		index = index - posLastGeneX
	return index - 1

##
## Main
##
if __name__ == '__main__':
		
	#parseOutputFile("test-orthoAlign.out")
	#parseOutputFile("test-duploss.out")
	#parseOutputFile("orthoAlign.out")
	#parseOutputFile("orthoAlign-transpo.out")
	
	#outputEvents("test-orthoAlign.out", "testevents.out")
	#outputEvents("orthoAlign-transpo.out", "testevents.out")
	#outputEvents("orthoAlign.out", "testevents.out")
	outputEvents("orthoAlign-bigTest.out", "testevents.out")
	#outputEvents("duploss-bigTest.out", "testevents.out")
	#outputEvents("duploss-bigTestMOD.out", "testevents.out")
	#outputEvents("orthoAlign-buggy.out", "testevents.out")