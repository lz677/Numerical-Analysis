# C:\Users\93715\python
# encoding: utf-8
"""
@author: LiuZhe
@license: (C) Copyright SJTU ME
@contact: LiuZhe_54677@sjtu.edu.cn
@file: SolveEquationSystem.py
@time: 2020/4/14 15:48
@desc:
"""

import numpy as np

# define A, b, x and Ax=b
n = 10
A = np.zeros([n, n])
b = np.zeros([n, 1])
x = np.ones([n, 1])
for i in range(n):
    if i == 0:
        A[i, 0] = 6
        A[i, 1] = 1
        b[i] = 7
    elif i == n - 1:
        A[i, -1] = 6
        A[i, -2] = 8
        b[i] = 14
    else:
        A[i, i - 1] = 8
        A[i, i] = 6
        A[i, i + 1] = 1
        b[i] = 15


# print(A)
# print(b)
# print(x)

def gauss_method(a_A, a_b, a_x):
    if a_A[0,0]