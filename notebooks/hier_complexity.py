# Author: Jason P. Smith
# Date: 13/10/22

import numpy as np
import math
import networkx as nx
import sys
import os


# Computes hierarchical complexity
# Input:  adjacency matrix of directed graph as numpy array
#         d1 and d2 specify whether to consider "in" "out" or "total" degrees
#         variability: takes arguments "std" (standard deviation), "var" (variance), or "alt" (which take the mean of the row and then std)
#         normalise: decised what to divide C[k] by, "none" (divides by 1), "k" (divides by k), "sqrt" (divides by sqrt(k)),
#                                                    "density" (divides by (1-density)*edges), "density2" (1-density)*edges/n
#         popstd: if 0 std divides by n, if 1 std divides by n-1
#         undirected: if True the inputed graph is considered as undirected, requires d1=d2='in', only used for normalise='density'
#         print_output: if True prints progress
# Output: (R,C) where R is the overall complextiy (float) and C is the vector of complexities for each degree value (vector of floats)
def hier_complexity(M, d1='in', d2='in', variability='std', normalise='density', popstd=0, undirected=True, print_output=True):
    if not print_output:
        sys.stdout = open(os.devnull, 'w')

    M = M.astype('uint8')
    print("Constructing Matrices.")
    n = len(M)

    density_factor = 1
    if undirected:
        M = M+np.transpose(M)
        M[M>0]=1
        print("Computing Degrees")
        degree = np.sum(M,axis=0)
        m = int(max(degree))
        print("Computing Degree Sequences")
        degree_seqs = [sorted([degree[i] for i in np.where(M[:,j])[0]]) for j in range(n)]
        print("Computing Degree Sets") #D[k] is the set of vertices of degree k, i.e D_k in the formula
        D = [np.where(degree==i)[0] for i in range(m+1)]
        num_edges = sum(degree)/2
        density = 2*num_edges/(n*(n-1))
    else:
        matrices = {'in':M, 'out':M.transpose()}
        if d1=='total' or d2=='total':
            matrices['total']=np.bitwise_or(matrices['in'],matrices['out'])
        print("Computing Degrees")
        degree1 = np.sum(matrices[d1],axis=0)
        if d1 == d2:
            degree2 = degree1
        else:
            degree2 = np.sum(matrices[d2],axis=0)
        m = int(max(degree1))
        print("Computing Degree Sequences")
        degree_seqs = [sorted([degree2[i] for i in np.where(matrices[d1][:,j])[0]]) for j in range(n)]
        print("Computing Degree Sets")#D[k] is the set of vertices of degree k, i.e D_k in the formula
        D = [np.where(degree1==i)[0] for i in range(m+1)]
        num_edges = sum(degree1)
        density = num_edges/(n*(n-1))

    C = [np.nan for i in range(m+1)]
    d = 0

    print(density)

    for k in range(1,m+1):                     #First summand of formula over D_k
        print("Considering D_"+str(k))
        if len(D[k])>1:
            d += 1
            # Compute matrix of sorted degree sequences
            A = np.array([degree_seqs[i] for i in D[k]])

            # Compute variability of columns
            if variability == 'std':
                V = np.std(A,axis=0,ddof=popstd)
            elif variability == 'var':
                V = np.var(A,axis=0)
            elif variability == 'alt':
                V = np.std(np.mean(A,axis=1),ddof=popstd)
            else:
                print("Error: invalid divideby")
                return -1

            # Divide variability by desired normalisation
            if normalise == 'none':
                C[k] = sum(V)
            elif normalise == 'k':
                C[k] = np.mean(V)
            elif normalise == 'sqrt':
                C[k] = sum(V)/math.sqrt(k)
            elif normalise == 'density':
                C[k] = sum(V)/((1-density)*num_edges)
            elif normalise == 'density2':
                C[k] = sum(V)/((1-density)*(num_edges/n))
                print(C[k])
            elif normalise == 'density3':
                C[k] = sum(V)/((1-density)*(num_edges/n)**2)
            elif normalise == 'density2std':
                C[k] = sum(V)/math.sqrt((1-density)*(num_edges/n))
            elif normalise == 'density3std':
                C[k] = sum(V)/math.sqrt((1-density)*(num_edges/n)**2)
            else:
                print("Error: invalid normalise")
                return -1

    R = np.nanmean(C)
    print("Number of distinct degrees is "+str(d))
    print("Total Complexity is "+str(R))

    if not print_output:
        sys.stdout = sys.__stdout__

    return (R,C)

# Input: The adjacency matrix of the graph of which the config model is required,
#        plus the keyword arguments of heir_complexity
# Output: heir complexity of the config model
def undirected_config_model_HC(M, variability='std', normalise='none', print_output=True):
    M = M.astype('uint8')
    degrees = sum(M)
    G = nx.configuration_model(degrees)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return hier_complexity(nx.to_numpy_array(G), 'in', 'in', variability=variability,
                            normalise=normalise, undirected=True, print_output=print_output)

# Input: The adjacency matrix of the graph of which the config model is required,
#        plus the keyword arguments of heir_complexity
# Output: heir complexity of the config model
def directed_config_model_HC(M, d1='in', d2='in', variability='std', normalise='none', print_output=True):
    M = M.astype('uint8')
    in_degree = sum(M)
    out_degree = sum(M.transpose())
    G = nx.directed_configuration_model(in_degree, out_degree)
    G = nx.DiGraph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    return hier_complexity(nx.to_numpy_array(G), d1, d2, variability=variability,
                            normalise=normalise, undirected=False, print_output=print_output)
