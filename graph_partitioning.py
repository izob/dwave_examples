
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import dwave.inspector
from dwave_qbsolv import QBSolv
from dwave.system.samplers import DWaveSampler           # Library to interact with the QPU
from dwave.system.composites import EmbeddingComposite 



# ! ------------- PARAMETERS D-WAVE -------------
sapi_token = 'DEV-5dc2cb00522e49c24b53362c103eeadc1ca8034c'
url = 'https://cloud.dwavesys.com/sapi/'


# ! ------------- PROBLEM PARAMETERS -------------
n = 6  # number of nodes
pr_edge = 0.3  # probability for edge creation
G = nx.gnp_random_graph(n=n, p=pr_edge, seed=123)  # returns a binomial graph
num_edges = len(G.edges)
delta = 0  # maximum number of edges that are incident to vertex
for i, degree in nx.degree(G):
    if delta < degree:
        delta =  degree
A = min([2*delta, n]) / 8

# * ------ QUBO
quadratic, linear = {}, {}
quadratic_dwave, linear_dwave = {}, {}
for i in range(n):
    for j in range(i+1, n):
        quadratic[i,j] = (2 * A)
        if (i,j) in G.edges:
            quadratic[i,j] += -1/2
        # quadratic_dwave[i,j] = 2*number_set[i]*number_set[j] / (2 * max1 * max2)
        quadratic[i,j] /=  (2 * A)
quadratic_dwave = quadratic.copy()

# # * ------ SOLVE PROBLEM (LOCAL)
response = QBSolv().sample_ising(linear, quadratic)
response_samples = list(response.samples())
response_energy = list(response.data_vectors['energy'])  # energies is the energy for each sample
response_counts = list(response.data_vectors['num_occurrences'])  # counts is the number of times each sample was found by qbsolv.
print('------------------ LOCAL ----------------------')
print(response_samples)
print(response_energy)


# * ------ SOLVE PROBLEM (DWAVE)
response_dwave = EmbeddingComposite(
    DWaveSampler(token=sapi_token, endpoint=url, solver='DW_2000Q_6')
    ).sample_ising(linear_dwave, quadratic_dwave, chain_strength=5, num_reads=10)
response_dwave_samples = list(response_dwave.samples())
response_dwave_energy = list(response_dwave.data_vectors['energy'])  # energies is the energy for each sample
response_dwave_counts = list(response_dwave.data_vectors['num_occurrences'])  # counts is the number of times each sample was found by qbsolv.
print('\n------------------ D-WAVE ----------------------')
print(response_dwave_samples)
print(response_dwave_energy)
print(response_dwave_counts)
dwave.inspector.show(response_dwave)

# nx.draw(G, pos=nx.circular_layout(G), with_labels=True)
# plt.show()


# import networkx as nx 
# from collections import defaultdict
# from itertools import combinations
# from dwave.system.samplers import DWaveSampler
# from dwave.system.composites import EmbeddingComposite
# import math

# # ------- Set tunable parameters -------
# num_reads = 1000
# gamma = 80

# # ------- Set up our QUBO dictionary -------

# # Initialize our Q matrix
# Q = defaultdict(int)

# # Fill in Q matrix
# for u, v in G.edges:
#     Q[(u,u)] += 1
#     Q[(v,v)] += 1
#     Q[(u,v)] += -2

# for i in G.nodes:
#     Q[(i,i)] += gamma*(1-len(G.nodes))

# for i, j in combinations(G.nodes, 2):
# 	Q[(i,j)] += 2*gamma

# print(Q)

# # ------- Run our QUBO on the QPU -------

# # Set chain strength
# chain_strength = gamma*len(G.nodes)


# response = QBSolv().sample_qubo(Q, chain_strength=chain_strength, num_reads=num_reads)

# # See if the best solution found is feasible, and if so print the number of cut edges.
# sample = response.record.sample[0]
# print(sample)

# # In the case when n is odd, the set may have one more or one fewer nodes
# if sum(sample) in [math.floor(len(G.nodes)/2), math.ceil(len(G.nodes)/2)]:
#     num_cut_edges = 0
#     for u, v in G.edges:
#         num_cut_edges += sample[u] + sample[v] - 2*sample[u]*sample[v]
#     print("Valid partition found with", num_cut_edges, "cut edges.")
# else:
#     print("Invalid partition.")