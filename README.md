# Error-decoding-community-detection

## Required Python Packages

1. itertools
2. numpy
3. networkx
4. random
5. copy
6. sys

## Running the script

The script requires three inputs

1. An edgelist representation of a graph. Each edge is separated by a space and on its own line
2. Number of times, x, to run the algorithm. Outputs x number of results.
3. Starting initial conditions. Must be 'Random' or 'Regular'. 


## Sample Command

```shell
nohup python -u Error_Decoding_Detection.py sample_graph.txt 5 Random &
```

## Output

1

('0', '2', '3', '4', '5', '6', '7', '1', '8', '9')

('13', '18', '15', '11', '17', '19', '14', '10', '12', '16')

('20', '28', '21', '25', '27', '22', '26', '23', '24', '29')

2

('0', '2', '3', '4', '5', '6', '7', '1', '8', '9')

('20', '28', '21', '25', '27', '22', '26', '23', '24', '29')

('13', '18', '15', '11', '17', '19', '14', '10', '12', '16')

3

('0', '2', '3', '4', '5', '6', '7', '1', '8', '9')

('20', '28', '21', '25', '27', '22', '26', '23', '24', '29')

('13', '18', '15', '11', '17', '19', '14', '10', '12', '16')

4

('20', '28', '21', '25', '27', '22', '26', '23', '24', '29')

('0', '2', '3', '4', '5', '6', '7', '1', '8', '9')

('13', '18', '15', '11', '17', '19', '14', '10', '12', '16')

5

('20', '28', '21', '25', '27', '22', '26', '23', '24', '29')

('0', '2', '3', '4', '5', '6', '7', '1', '8', '9')

('13', '18', '15', '11', '17', '19', '14', '10', '12', '16')




