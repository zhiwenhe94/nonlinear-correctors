# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import sympy as sp
import itertools
import math
import random
from math import comb
from collections import defaultdict


def combina(n,k):
    return comb(n,k)

# For given parameters n and t, create matrix H
def generate_matrix_H(n, t):
    # Create a zero matrix of size (n+1) x (n+1)
    H = np.zeros((n + 1, n + 1), dtype=int)

    for s in range(0,n+1):
        for d in range(0,t+1):
                H[d, s] = combina(s, d) * combina(n - s, t - d)

    return H

# For given parameters n, m, and t, create vector S1
def generate_vector_S1(n,m,t):
    # Calculate the constant factor
    factor=2**(n-m-t)*combina(n,t)

    # Create a zero vector of length (n+1)
    S1=np.zeros(n+1,dtype=int)

    # Fill the first t+1 elements
    for k in range(t+1):
        S1[k]=factor*combina(t,k)

    return S1

def generate_vector_S2(n):
    # Create a zero vector of length (n+1)
    S2=np.zeros(n+1,dtype=int)
    # Fill the first n+1 elements
    for k in range(n+1):
        S2[k]=combina(n,k)

    return S2


def generate_vectors(n):
    # Calculate the binomial coefficient C_n^s for each s
    binomials = [math.comb(n, s) for s in range(n + 1)]
    # For each s, generate all possible values in the range [0, C_n^s]
    ranges = [range(binomials[s] + 1) for s in range(n + 1)]
    # Generate all possible vectors
    all_vectors = [np.array(vec) for vec in itertools.product(*ranges)]

    return all_vectors

def generate_solutions(n,m,t):
    solutions_set=[]
    # Generate the coefficient matrix H
    H=generate_matrix_H(n, t)
    S1=generate_vector_S1(n,m,t)
    vectors=generate_vectors(n)
    for v in vectors:
        #if H*v==S1
        if np.array_equal(np.dot(H,v),S1):
            solutions_set.append(v.tolist())

    return solutions_set

# Select 2**(n-m) solutions from the solutions_set with replacement, and sum them to equal S2
# Generate all possible selections from the solutions_set and filter out those that sum to S2
def find_all_solutions_sum_s2(solutions_set, n, m):
    target_count = 2 ** (m)
    S2 = generate_vector_S2(n)

    # Generate all possible combinations
    all_combinations = itertools.combinations_with_replacement(solutions_set, target_count)

    valid_combinations = []
    for combo in all_combinations:
        sum_vector = np.sum(combo, axis=0)
        if np.array_equal(sum_vector, S2):
            valid_combinations.append(combo)

    return valid_combinations



# Generate the n-dimensional vector space over F_2
def generate_f2_vector_space(n):
    # Create the finite field F_2
    field=sp.GF(2)
    # Generate the n-dimensional vector space over F_2
    vectors=[]
    for i in range(2**n):
        # Convert integer i into a binary vector
        vector=[int(x) for x in bin(i)[2:].zfill(n)]
        vectors.append(sp.Matrix(vector))

    return vectors

#Define the vector weight
def vector_weight(vector):
    return sum(x==1 for x in vector)

# Partition the n-dimensional vector space over F_2 based on weight
def group_vectors_by_weight(vectors):
    weight_groups=defaultdict(list)
    for vector in vectors:
        weight=vector_weight(vector)
        weight_groups[weight].append(vector)
    return weight_groups

# According to the combinations generated in find_all_solutions_sum_s2,
# partition the n-dimensional vector space over F_2 into 2**m parts,
# where the size of each part is 2**(n-m) and the number of vectors with different weights in each part corresponds to the vectors in the combination generated in find_all_solutions_sum_s2.


def random_partition_by_weight(weight_groups,combo,n):
    Partition=[]
    group1_sizes=np.zeros(n+1,dtype=int)
    for group_sizes in combo:
        group=[]
        for weight,vectors in weight_groups.items():
            random.shuffle(vectors)
            if group1_sizes[weight]+group_sizes[weight] > len(vectors):
                raise ValueError("Not enough vectors to satisfy the partition for this weight")
            group.extend(vectors[group1_sizes[weight]:group_sizes[weight]+group1_sizes[weight]])
            group1_sizes[weight]+=group_sizes[weight]
        # Convert each matrix in the group to a tuple and store it in Partition
        group_as_tuples={tuple(matrix) for matrix in group}
        Partition.append(group_as_tuples)

    return Partition

# examples
n = 6
m = 3
t = 2

solutions_set = generate_solutions(n, m, t)
all_valid_combinations = find_all_solutions_sum_s2(solutions_set, n, m)
print(f"Number of valid combinations: {len(all_valid_combinations)}")
for idx, combo in enumerate(all_valid_combinations):
    print(f"Combination {idx + 1}:")
    for sol in combo:
        print(sol)


vectors=generate_f2_vector_space(n)
weight_groups=group_vectors_by_weight(vectors)
print(f"Number of valid combinations: {len(all_valid_combinations)}")
for idx, combo in enumerate(all_valid_combinations):
    print(f"The {idx + 1}th solution of F(x):")
    Partition=random_partition_by_weight(weight_groups,combo,n)
    print(Partition)





# See PyCharm help at https://www.jetbrains.com/help/pycharm/
