
######################################
###Self Global Alignment Functions####
######################################

######################################################
# findOrthologsBySelfGlobalAlignment
# Parameters:
# Description: Performs a self global alignment to identify duplicate operons
######################################################
def findOrthologsBySelfGlobalAlignment(strain, coverageTracker):
    print('Performing self-global alignment on strain: %s' %(strain.name))
    lossEvents = []
    duplicationEvents = []
    fragments = strain.genomeFragments
    
    for x in range(0, coverageTracker):
        if coverageTracker[x] == False:
            bestScore = 1000    #Make the best score some large numer
            bestEvent = None    #Initialize the event
            minDistance = 1000  #Used to track the minimum distance from singleton to operon that has an identical gene
            
            filteredList = iter(filter(lambda x: x.fragmentIndex == x, fragments)) #Get the fragment we need based on the index
            foundFragment = next(filteredList, None)
            
            if foundFragment == None:
                print('Something went wrong! Unable to find specific genome fragment!')
            else:
                if len(foundFragment.sequence) > 1: #We're processing an operon
                    print('')
                    
                else: #We're processing a singleton
                    print('')
    return None