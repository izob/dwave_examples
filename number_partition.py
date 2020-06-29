import random
import numpy as np
import dimod
import dwave.inspector
from dwave_qbsolv import QBSolv
from dwave.system.samplers import DWaveSampler           # Library to interact with the QPU
from dwave.system.composites import EmbeddingComposite 



# ! ------------- PARAMETERS D-WAVE -------------
sapi_token = 'DEV-5dc2cb00522e49c24b53362c103eeadc1ca8034c'
url = 'https://cloud.dwavesys.com/sapi/'


# ! ------------- PROBLEM PARAMETERS -------------
# np.random.seed(seed=1234567890)
# number_set = np.random.choice(range(5),10, replace=True)
# number_set[3] = 2
number_set = np.array([1, 2, 3])
max1, max2 = np.sort(number_set)[-1], np.sort(number_set)[-2]

# * ------ QUBO
quadratic, linear = {}, {}
quadratic_dwave, linear_dwave = {}, {}
linear_dwave[0] = 1
linear[0] = 0
for i in range(len(number_set)):
    for j in range(i+1, len(number_set)):
        quadratic[i,j] = 2*number_set[i]*number_set[j] 
        quadratic_dwave[i,j] = 2*number_set[i]*number_set[j] / (2 * max1 * max2)

# * ------ SOLVE PROBLEM (LOCAL)
response = QBSolv().sample_ising(linear, quadratic)
response_samples = list(response.samples())
response_energy = list(response.data_vectors['energy'])  # energies is the energy for each sample
response_counts = list(response.data_vectors['num_occurrences'])  # counts is the number of times each sample was found by qbsolv.
print('------------------ LOCAL ----------------------')
print(number_set)
print(response_samples)
print(response_energy)

# * ------ SOLVE PROBLEM (DWAVE)
response_dwave = EmbeddingComposite(
    DWaveSampler(token=sapi_token, endpoint=url, solver='DW_2000Q_6')
    ).sample_ising(linear_dwave, quadratic_dwave, chain_strength=100, num_reads=10)
response_dwave_samples = list(response_dwave.samples())
response_dwave_energy = list(response_dwave.data_vectors['energy'])  # energies is the energy for each sample
response_dwave_counts = list(response_dwave.data_vectors['num_occurrences'])  # counts is the number of times each sample was found by qbsolv.
print('\n------------------ D-WAVE ----------------------')
print(number_set)
print(response_dwave_samples)
print(response_dwave_energy)
print(response_dwave_counts)
dwave.inspector.show(response_dwave)
