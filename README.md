# Error-decoding-community-detection

## Required Python Packages

1. itertools
2. numpy
3. networkx
4. random
5. copy
6. sys
7. csv

## Running the script

The script requires three inputs

1. An edgelist representation of a graph. Each edge is separated by a space and on its own line
2. Number of times, x, to run the algorithm. Outputs x number of results.
3. Starting initial conditions. Must be 'Random' or 'Regular'. 


## Sample Command

nohup python -u Error_Decoding_Detection.py edgelist.txt 10 Random

