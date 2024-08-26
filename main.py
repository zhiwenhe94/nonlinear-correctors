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

#对于给定的参数n，t创建矩阵H
def generate_matrix_H(n, t):
    # 创建一个大小为 (n+1) x (n+1) 的零矩阵
    H = np.zeros((n + 1, n + 1), dtype=int)

    for s in range(0,n+1):
        for d in range(0,t+1):
                H[d, s] = combina(s, d) * combina(n - s, t - d)

    return H

#对于给定参数n,m,t创建向量S1
def generate_vector_S1(n,m,t):
    #计算常数因子
    factor=2**(n-m-t)*combina(n,t)

    #创建一个长度为(n+1)的零向量
    S1=np.zeros(n+1,dtype=int)

    #填充前t+1个元素
    for k in range(t+1):
        S1[k]=factor*combina(t,k)

    return S1

def generate_vector_S2(n):
    #创建一个长度为(n+1)的零向量
    S2=np.zeros(n+1,dtype=int)
    #填充前n+1个元素
    for k in range(n+1):
        S2[k]=combina(n,k)

    return S2


def generate_vectors(n):
    # 计算每个s对应的C_n^s组合数
    binomials = [math.comb(n, s) for s in range(n + 1)]
    # 为每个s生成介于[0, C_n^s]的所有可能值
    ranges = [range(binomials[s] + 1) for s in range(n + 1)]
    # 生成所有可能的向量
    all_vectors = [np.array(vec) for vec in itertools.product(*ranges)]

    return all_vectors

def generate_solutions(n,m,t):
    solutions_set=[]
    #生成系数矩阵H
    H=generate_matrix_H(n, t)
    S1=generate_vector_S1(n,m,t)
    vectors=generate_vectors(n)
    for v in vectors:
        #如果H*v==S1
        if np.array_equal(np.dot(H,v),S1):
            solutions_set.append(v.tolist())

    return solutions_set

#从solutions_set中有重复地挑选2**(n-m)个解，相加等于S2
# 从solutions_set中生成所有可能的挑选方法，并过滤出和等于S2的组合
def find_all_solutions_sum_s2(solutions_set, n, m):
    target_count = 2 ** (m)
    S2 = generate_vector_S2(n)

    # 生成所有可能的组合
    all_combinations = itertools.combinations_with_replacement(solutions_set, target_count)

    valid_combinations = []
    for combo in all_combinations:
        sum_vector = np.sum(combo, axis=0)
        if np.array_equal(sum_vector, S2):
            valid_combinations.append(combo)

    return valid_combinations



#生成F_2上的n维向量空间
def generate_f2_vector_space(n):
    #创建有限域F_2
    field=sp.GF(2)
    #生成F_2上的n维向量空间
    vectors=[]
    for i in range(2**n):
        #将整数i转换维二进制形式的向量
        vector=[int(x) for x in bin(i)[2:].zfill(n)]
        vectors.append(sp.Matrix(vector))

    return vectors

#定义向量weight
def vector_weight(vector):
    return sum(x==1 for x in vector)

#将F_2上的n维向量空间按照weight进行划分
def group_vectors_by_weight(vectors):
    weight_groups=defaultdict(list)
    for vector in vectors:
        weight=vector_weight(vector)
        weight_groups[weight].append(vector)
    return weight_groups

#按照find_all_solutions_sum_s2中生成的combiantion，
# 将F_2上的n维向量空间进行划分为2**m个parts，
# 每个part的大小为2**(n-m)且每个part中的不同weight的向量个数对应find_all_solutions_sum_s2中生成的combiantion中的向量

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
        #将group中每个矩阵转换维tuple并存储在Partition中
        group_as_tuples={tuple(matrix) for matrix in group}
        Partition.append(group_as_tuples)

    return Partition

# 示例代码
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
    print(f"F(x)的第 {idx + 1}种解:")
    Partition=random_partition_by_weight(weight_groups,combo,n)
    print(Partition)



#https://github.com/zhiwenhe94/nonlinear-correctors.git


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
