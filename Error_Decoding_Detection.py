"""
Linear algorithm for community Detection

Error-Correcting Decoders for Communities in Networks

Krishna Bathina and Filippo Radicchi
"""

# imports used in function
import itertools
import numpy as np
import networkx as nx
import random
import copy
import sys


def F(a,x): 
	"""
	Helper function that takes two inputs. 

	a is the multiplicative factor
	x is the input for the hyperbolic tangent
	"""

	atanh = a*np.tanh(x/2.0)
	# if atanh = 1, return a large value
	if atanh == 1.0: return np.log(100)

	# if atanh = -1, return a small value
	elif atanh == -1.0: return np.log(0.0001)

	# else, return the true value
	else: return np.log((1.0 + atanh)/(1.0 - atanh))

def parity_check(graph, xi_f, Z, edges, l_node, a, b, Ns,P_ratio,  time = 'continue'): #################
	"""
	Function that makes solves for the parity bit and information bit variables.

	graph = networkx graph
	xi_f = dictionary of xi for each edge
	Z = constant value of Z
	edges = list of edges
	l_node = llr of nodes
	a = constant used in F(a,x)
	b = constant used in F(a,x)
	Ns = constant used to in parity check
	P_ratio = ratio of P_in and P_out
	time = iteration number
	"""

	# for the nodes, calculate the final LLR
	final_node = l_node
	for i in l_node:
		final_node[i] += (Ns[i]+1)*F(b,Z) + sum(F(a,xi_f[(s,i)]) for s in graph.neighbors(i))
	
	# if LLR is less than 0: the community for the node is 1 
	# if LLR is greater than 0: the community for the node is 0
	final_node = {i:1 if final_node[i] < 0 else 0 for i in final_node}
	final_edge = {(i,j):0 for i,j in edges}
	
	# update edges in the same way
	for i,j in edges:
		# precalculate the multiplicative factor
		tanh = np.tanh(xi_f[(i,j)]/2) * np.tanh(xi_f[(j,i)]/2)
		try:
			# try calculating the final edge LLR
			final_edge[(i,j)] = np.log(P_ratio) + np.log((1.0+tanh)/(1.0-tanh))
		except ZeroDivisionError:
			# if denominator is 0, set LLR for the edge to be large
			final_edge[(i,j)] = 100

	# similar logic to above
	# if LLR is less than 0: pairty bit for the edge is 1 
	# if LLR is greater than 0: parity bit for the edge is 0
	final_edge = {edge:1 if final_edge[edge] < 0 else 0 for edge in edges}
	
	if time == 'continue': # used to check if parity equations are solved
		for i,j in final_edge: # do parity check on each edge
			if (final_node[i] + final_node[j] + final_edge[(i,j)]) % 2 != 0:
				return False # return false if any fail, algorithm keeps running
		# if all edges pass parity check, return the communities and parity bits
		return final_node, final_edge

	# in the final check, return the communities and parity bits
	elif time == 'final': return final_node, final_edge
	
def BP(graph, initial_condition): 
	"""
	function that calculates LLRS using BP algorithm

	graph = networkx graph
	initial_condition = starting conditions
	"""
	# if there are no edges, return remaining nodes in one community
	if len(graph.edges()) == 0: return [{x:0 for x in graph.nodes()},()]

	nodes = list(graph.nodes()) # list of nodes
	N = len(graph) # number of nodes in the graph

	# calculation for average degree: sum of total number of neighbors divided by the number of nodes
	k = sum(dict(graph.degree(graph.nodes())).values()) / len(graph)

	# set constant value to ensure that Pin is greater than Pout
	alpha = 1.2

	# solve for Pin and Pout
	P_in = alpha * (k + k**0.5 ) / N
	P_out = max(0, 2*k/N - P_in )
	# if P_out is zero, set to small value
	if P_out == 0: P_out = 0.1

	# get list of all edges (back and forth)    
	edges = list(graph.edges()) # number of edges
	edges.extend([q[::-1] for q in edges]) # number of directed edges
	edges = set(edges)
	
	# total number of edges that don't exist
	num_non_edges = N * (N-1) - len(edges)
	
	# constants used in f(a,x). Calculating now for use later
	a = (P_in - P_out)/(P_in + P_out)
	b = (P_out - P_in) / (2.0 - P_in - P_out)

	# if starting initial conditions are random
	if initial_condition == 'Random':
		# make an LLR dictionary with values from a uniform distribution from -1 to 1
		l_node = {x:random.uniform(-1,1) for x in nodes} 

	# if starting initial conditions are regular
	elif nitial_condition == 'Regular':
		# make an LLR dictionary setting all LLR to zero
		l_node = {x:0.0 for x in nodes} 
		# set a random node's LLR to be 1
		l_node[random.choice(nodes)] = 1

	# initializing values of xi
	xi_i = {(i,j):l_node[i] for i,j in edges}

	# constants used for later
	## numerator used for Z
	numerator = sum((N-1-graph.degree(x)) * l_node[x] for x in nodes)
	## constant similar to numberator with the LLR
	Ns = {x:N-1-graph.degree(x) for x in nodes}

	# Calculating Z
	## If there are non-edges that exist - calcualte Z normally
	if num_non_edges: Z = numerator/num_non_edges # initial value of z
	## else, assume small value of denominator 
	else: Z = 10 * numerator

	# calculating ratio of P_in and P_out
	p_ratio = P_in/P_out


	t = 0 # intiailize time
	error = 1 # initialize error

	# run while loop for less than 1000 iterations or while the error is large
	while error > 10**-14 and t < 1000:

		# update xi
		## first account for the LLR
		xi_f = {(i,j):l_node[i] for i,j in edges}
		## update for each edge
		for i,j in edges:
			xi_f[(i,j)] += Ns[i]*F(b,Z) + sum(F(a,xi_i[(s,i)]) for s in graph.neighbors(i) if s!= j)
			
		# for every 10 iterations
		if not t%10: 
			# return results of parity check - continuous
			finals = parity_check(graph, xi_f, Z, edges, l_node, a, b, Ns, p_ratio)
			# if all parity check equations pass, return the 
			if finals: return finals  
			
		# save value of xi used for next iteration
		xi_i = copy.deepcopy(xi_f)
		# save value of Z used for next iteration
		Z_old = Z

		# if there are non-edges 
		if num_non_edges: 
			# update Z 
			Z = F(b,Z) + (numerator + sum(F(a,xi_i[(i,j)]) for i,j in edges)) / num_non_edges
		# if there are no non-edges
		else: 
			# make denominator very small and calculate Z
			Z = F(b,Z) + 10 * (numerator + sum(F(a,xi_i[(i,j)]) for i,j in edges))
		# calculate error as the difference of Z
		error = abs(Z-Z_old)
			
		# increment time
		t += 1  

	# return results of parity check - final time         
	return parity_check(graph, xi_f, Z, edges, l_node, a, b, Ns, p_ratio, time = 'final')

def subgraph_run(G, nodes, initial_condition):
	"""
	function that runs sets up and calls the BP algorithm for subgraphs of G that can be further split

	G = original networkx graph
	nodes = list of nodes that make of the subgraph of the original graph G
	initial_condition = starting conditions
	"""

	# make an empty network, add nodes from the nodelist, and add edges from the original graph
	T=nx.Graph() 
	T.add_nodes_from(nodes)
	T.add_edges_from([x for x in list(nx.edges(G, nodes)) if x[0] in nodes and x[1] in nodes])

	# run BP on the sub graph and save the values
	final_node, final_edge = BP(T, initial_condition)
	one = tuple(x for x in final_node if final_node[x] == 0) # make a tuple of nodes with a community membreship of 0
	two = tuple(x for x in final_node if final_node[x] == 1) # make a tuple of nodes with a community membreship of 0

	return one,two # return the communities

def run(G,initial_condition):  
	"""
	function that sets up and calles the BP algorithm

	G = networkx object, undirected graph
	initial_condition = starting conditions
	"""

	final_node = {} # dictionray that will hold final community values of nodes (0 or 1)
	t = 0 # time increment variable

	# while loop that only stops if there is a split (2 communities)
	# if a split doesn't occur after 10 tries, set the final node and final edge values

	while 1 not in final_node.values() or 0 not in final_node.values():

		if t == 10: # break after 10 attempts
			break 

		final_node, final_edge = BP(G, initial_condition) # run BP and save values
		t += 1 # increment the time variable

	first = tuple(x for x in final_node if final_node[x] == 0) # make a tuple of nodes with a community membreship of 0
	second = tuple(x for x in final_node if final_node[x] == 1) # make a tuple of nodes with a community membreship of 1

	 # if all ten attempts were used, return a list of nodes separated by community. Because there was no split, 
	 # one list is empty and the other contains all the nodes
	if t == 10: 
		return [first,second]
	
	# if t != 10
	else:
		# make an initial results dictionary. The keys in the dictionary are lists of nodes in the community and the 
		# value indicates if the community can be further split (1 = yes, 0 = no)
		results = {first:1,second:1} # initial values indicate that the communities can possibly be further split
		
		while 1 in results.values(): # run until no more communities can be further split

			previous_results = results # store previous results

			# update communities that cannot be split into results dictionary
			results = {x:0 for x in previous_results if not previous_results[x]}

			# for each community that can be split
			for com in (x for x in previous_results if previous_results[x]):
				a,b = subgraph_run(G,com, initial_condition) # attempt to split the community by running BP on the smaller network
				# if any of the splits are empty, the original com is done, set to 0s
				if not a or not b: results[com] = 0
				# communities can keep on splitting
				else: results[a] = results[b] = 1
					
		return list(results.keys()) # return list of lists of communities


if __name__ == "__main__":
	
	edgelist = sys.argv[1] # edgelist
	iterations = int(sys.argv[2]) # number of iterations (will output results from this many runs of the function)
	initial_condition = sys.argv[3] # number of iterations (will output results from this many runs of the function)

	# exit if the proper initial conditions are not inputted
	if initial_condition != 'Random' and initial_condition != 'Regular':
		sys.exit("Incorrect initial conditions. Must be 'Random' or 'Regular'")

	# exit if not enough initial conditions (must be greater than 0)
	if iterations < 1: 
		sys.exit("Number of runs must be greater than 0.")

	# exit if graph_files are not properly formatted 
	try:
		# read edgelist using networkx 
		graph = nx.read_edgelist(edgelist)
	except:
		sys.exit("Edgelist not properly configured.") 

	# place to store results
	results = []
	for i in range(iterations): # for each iteration
		print('iteration =', i) # print iteration number
		results.append(run(graph,initial_condition)) # store community labels in results

	with open("results.csv","w") as f: # open save file
		for x in range(len(results)):
			f.write(str(x+1) + '\n')
			for com in results[x]:
				f.write(str(com) + '\n')