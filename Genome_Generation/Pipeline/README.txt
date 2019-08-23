Generator : Generates simulated gene strains within an ancestral tree by creating a root ancestor genome and producing events that affect the ancestral genome down different paths of the tree. The output is a separate text file containing a gene order for every leaf node within the tree. The output also includes a text file for the gene strain of the root ancestor.

Author : Meghan Chua
Contact : umchua2@myumanitoba.ca

Bioinformatics Laboratory
Computer Science Department
University of Manitoba

USAGE NOTES :
-------------------------------------------------------------------

python genomeGenerator.py  tree length operons events -dliste

Parameters description :

tree is the file containing the tree structure we are trying to generate data for (must be newick format)

length is the number of genes we would like for the root ancestor

operons is the minimum number of operons we want the root ancestor to have

events is the number of events we want to occur per branch in the tree

-d allows duplication events to occur, two parameters (probability, p-value).

-l allows loss events to occur, two parameters (probability, p-value).

-i allows inversion events to occur, two parameters (probability, p-value).

-s allows substitution events to occur, one parameters (probability).

-t allows transposition events to occur, two parameters (probability, p-value).

-e special case that sets equal distribution of allowed events per branch, and only allows one inversion to occur in the whole tree. Events must be a multiple of num-allowed-events (eg. events = 8 with -d 0.20 0.7 -l 0.20 0.7 -i 0.20 0.7 -t 0.20 0.7, num-allowed-events = 4, each event will occur 2 times per branch). Will always include duplications, losses, and inversions, but substitutions are optional.

*probability (0.0 - 1.0) : Probability of that event to occur. The total probability of all events we allow must equal 1.0.

**p-value (0.0 - 1.0) : Parameter that influences the size of the events (number of operons/genes that will be affected). The larger the p-value, the smaller the event size.

EXAMPLE :
------------------------------------------------------------------

python genomeGenerator.py tree2Leaf.dnd 120 5 4 -d 0.20 0.7 -l 0.20 0.7 -i 0.20 0.7 -t 0.20 0.7 -s 0.20 -e

Input and output for this example can be found in the zip folder exampleTest.zip