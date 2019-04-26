import copy
import globals

######################################################
# computeOperonArrangements
# Parameters:
# Description: Takes a list of events and computes Forward Conserved, Transposed Forward, Inverted, Inverted Transposed, and Lost regions
######################################################
def computeOperonArrangements(events):
    
    #Various regions of the genome
    lostRegions = []
    
    #Sort events by x coordinates
    eventsCopy = copy.deepcopy(events)
    eventsCopy.sort(key=lambda x:x.fragmentDetails1.fragmentIndex, reverse=False)
    
    while len(eventsCopy) > 0:
        currEvent = eventsCopy.pop(0)
        
        if currEvent.score == - 1:          #This is an operon that has been lost
            lostRegions.append(currEvent)
        else:                               #This an operon on the dot plot
            consecutiveRegion = [] #Stores a consecutive region
            consecutiveRegion.append(currEvent)           
            foundNeighbor = True #Keeps track of whether we found a neighboring operon
            
            while foundNeighbor: #Continue looping as long as we keep on finding a neighbor
                foundNeighbor = False #Reset the tracker
                
                prevEvent = consecutiveRegion[len(consecutiveRegion)-1] #Gets the last operon in the the consecutive operons list
                currEvent = eventsCopy[0] #Get the next available operon
                
                yDistance = abs(currEvent.fragmentDetails2.fragmentIndex - prevEvent.fragmentDetails2.fragmentIndex) #The distance on the y-axis
                
                if yDistance < globals.yDistanceThreshold: #If the y-Distance is less than the threshold then add it to the consecutive region
                    consecutiveRegion.append(eventsCopy.pop(0))
                    foundNeighbor = True
                
            
    
    
    return None