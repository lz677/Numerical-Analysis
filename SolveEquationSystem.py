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
# n = 10
# n = 30
n = 100
A = np.zeros([n, n])
b = np.zeros([n, 1])
x = np.ones([n, 1])

for matrix_i in range(n):
    if matrix_i == 0:
        A[matrix_i, 0] = 6
        A[matrix_i, 1] = 1
        b[matrix_i] = 7
    elif matrix_i == n - 1:
        A[matrix_i, -1] = 6
        A[matrix_i, -2] = 8
        b[matrix_i] = 14
    else:
        A[matrix_i, matrix_i - 1] = 8
        A[matrix_i, matrix_i] = 6
        A[matrix_i, matrix_i + 1] = 1
        b[matrix_i] = 15

print("det(A) = " + str(np.linalg.det(A)))

# print(A)
# print(b)
# print(x)

"""
高斯顺序消元
"""


def gauss_method(a_gauss, b_gauss, x_gauss):
    for k in range(n):
        for i in range(k + 1, n):
            if a_gauss[k, k]:
                m = a_gauss[i, k] / a_gauss[k, k]
                a_gauss[i, k] = 0
                b_gauss[i] = b_gauss[i] - m * b_gauss[k]
                for j in range(k + 1, n):
                    a_gauss[i, j] = a_gauss[i, j] - m * a_gauss[k, j]
            else:
                print("a" + str(k) + str(k) + "=" + str(a_gauss[k, k]))
    for i in range(n - 1, -1, -1):
        sum = 0
        if i == n - 1:
            x_gauss[i] = b_gauss[i] / a_gauss[i, i]
        else:
            for j in range(i + 1, n):
                sum = sum + a_gauss[i, j] * x_gauss[j]
            x_gauss[i] = (b_gauss[i] - sum) / a_gauss[i, i]


# gauss_A = A.copy()
# gauss_b = b.copy()
# x_gauss_app = np.zeros([n, 1])
# gauss_method(gauss_A, gauss_b, x_gauss_app)
# print(x-x_gauss_app)
# np.savetxt('delta_x.txt', x-x_gauss_app, fmt='%0.16f')
# print(x - x_gauss_app)

"""
追赶法
Ax = f;
A = |b1 c1                | = |alpha1                    | |1 beta_1            |
    |a2 b2 c2             |   |r2     alpha2             | |  1  beta_2         |
    |0  ... ... ...       |   |         ...    ......    | |    ...             |
    |      an-1 bn-1 cn-1 |   |                          | |           beta_n-1 |
    |           an    bn  |   |              rn   alpha_n| |           1        |
公式：
1.求beta
    beta_1 =  c1/b1
    beta_i = ci/(bi-ai*beta_i-1) i = 2,3,...n-1
2.Ly =f 
    y1 = f1/b1
    yi = (fi-ai*yi-1)/(bi-ai*beta_i) i = 2,3,...n
3.Ux = y
    xn = yn
    xi = yi - beta_i* xi+1  i = n-1, n-2,....,2,1
    
"""


def chasing_method(a_chasing, f_chasing):
    beta = np.zeros([n - 1, 1])
    # 计算bata_i
    # beta_1 = c1/b1
    beta[0, 0] = a_chasing[0, 1] / a_chasing[0, 0]
    for i in range(n - 1):
        # beta_i = ci/(bi-ai*beta_i-1) i = 2,3,...n-1
        beta[i, 0] = a_chasing[i, i + 1] / (a_chasing[i, i] - a_chasing[i, i - 1] * beta[i - 1, 0])

    # 解方程Ly=f  用y冲掉f
    # f1 <-- y1 = f1/b1
    f_chasing[0, 0] = f_chasing[0, 0] / a_chasing[0, 0]
    for i in range(1, n):
        # fi <-- yi = (fi-ai*yi-1)/(bi-ai*beta_i) i = 2,3,...n
        f_chasing[i, 0] = (f_chasing[i, 0] - a_chasing[i, i - 1] * f_chasing[i - 1, 0]) / (
                a_chasing[i, i] - a_chasing[i, i - 1] * beta[i - 1, 0])

    # 解方程 Ux=y 用x冲掉f
    # xn = yn
    # do nothing

    # xi = yi - beta_i* xi+1  i = n-1, n-2,....,2,1
    for i in range(n - 2, -1, -1):
        f_chasing[i, 0] = f_chasing[i, 0] - beta[i] * f_chasing[i + 1, 0]


chasing_a = A.copy()
chasing_f = b.copy()
chasing_method(chasing_a, chasing_f)
# print(chasing_f)
# np.savetxt('delta_x.txt', x - chasing_f, fmt='%0.32f')

# x = np.linalg.solve(A, b, )
"""
列主元素消去法
Ax=b.本算法用A的具有航交换的列主元素消去法，
消元结果冲掉A,乘数mij冲掉aij，计算x冲掉常数项b，行列式存放在det中
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
记录一次能气死的bug： 由于直接用b[k] 它是一个矩阵，在替换的时候出现了 
buf_b = b[k]
b[k] = b[k+1]
b[k+1] = buf_b
由于是矩阵的原因，三个是公用地址的关系，所以并没有起到交换b[k]与b[k+1]的作用
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!！！！！
"""


def main_column_method(a_m_col, b_m_col):
    det = 1
    for k in range(n - 1):
        max_value = 0
        ik_max_sub = 0
        # 按列选主元
        for i in range(k, n - 1):
            if i == k:
                max_value = abs(a_m_col[i, k])
                ik_max_sub = i
            # print("max_value in " + str(max_value))

            if abs(a_m_col[i + 1, k]) > abs(a_m_col[i, k]):
                max_value = abs(a_m_col[i + 1, k])
                # print("max_value bijiao " + str(max_value))
                ik_max_sub = i + 1

        # 如果最大列主元>0 停止计算 然后令行列式 = 0
        if max_value == 0:
            det = 0
            print('max_value is zero')
            return -1

        # 如果该主元不是最大主元 则交换行
        if k != ik_max_sub:
            # bk <--> bik
            buf_b = b_m_col[k, 0]
            b_m_col[k, 0] = b_m_col[ik_max_sub]
            b_m_col[ik_max_sub, 0] = buf_b
            # det = -det
            det *= -1
        # akj = aik_j (j= k, k+1, ... n)
        for j in range(k, n):
            buf = a_m_col[k, j]
            a_m_col[k, j] = a_m_col[ik_max_sub, j]
            a_m_col[ik_max_sub, j] = buf

        # 消元计算
        for i in range(k + 1, n):

            # aik <- mik = aik/akk
            a_m_col[i, k] = a_m_col[i, k] / a_m_col[k, k]

            # aij = aij - mik * akj  j=k+1,...,n
            for j in range(k + 1, n):
                a_m_col[i, j] = a_m_col[i, j] - a_m_col[i, k] * a_m_col[k, j]

            # bi = bi-mik*bk
            b_m_col[i, 0] = b_m_col[i, 0] - a_m_col[i, k] * b_m_col[k, 0]

        det *= a_m_col[k, k]

    # 如果ann = 0， 停止计算 det=0； 否则 det = det * ann;
    if a_m_col[n - 1, n - 1] == 0:
        det = 0
        print('ann is zero')
        return -1
    else:
        det *= a_m_col[n - 1, n - 1]
        print(det)

    # bn = bn/ann
    # bi = [bi - sum aij*bj(j=i+1,...n)]/aii  i = n-1, n-2 ,.....,1
    for m_col_back_i in range(n - 1, -1, -1):
        m_col_back_sum = 0
        if m_col_back_i == n - 1:

            b_m_col[m_col_back_i, 0] = b_m_col[m_col_back_i, 0] / a_m_col[m_col_back_i, m_col_back_i]
        else:
            for m_col_back_j in range(m_col_back_i + 1, n):
                m_col_back_sum += a_m_col[m_col_back_i, m_col_back_j] * b_m_col[
                    m_col_back_j]
            b_m_col[m_col_back_i] = (b_m_col[m_col_back_i] - m_col_back_sum) / a_m_col[
                m_col_back_i, m_col_back_i]


m_col_A = A.copy()
m_col_x = b.copy()
main_column_method(m_col_A, m_col_x)
# print(x-m_col_x)
# print(m_col_A)
np.savetxt('delta_x.txt', x - m_col_x, fmt='%0.90f')

"""
Jacobi 迭代
计算公式：
x0 = (x0_1, x0_2,.....x0_n)^T
xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数

终止标准：采用后验估计
||x* - xk||<= q/(q-1) ||xk - xk+1||
q = ||B|| = max_lamda_B

A = |a1_1                    |    |0                         |   |0 a1_2 ...  ...  a1_n-1  a1_n  |
    |    a2_2                |    |a2_1  0                   |   |  0  ...         a2_n-1  a2_n  |
    |         a3_3           | +  |  .      ...              | + |           ...                 |
    |              ...       |    |an-1_1          ...       |   |                   0     an-1_n|
    |                    an_n|    |an_1            an_n-1   0|   |                           0   |
A = D-L-U
B = I - inv(D)*A
"""


# p_significant_digits=-1时使用迭代次数作为终止标准，p_significant_digits！=-1 时采用后验估计和最大迭代次数（默认10000）做终止标准
def jacobi_iterative_method(a_jacobi, b_jacobi, x_jacobi, p_significant_digits=-1, max_iterative_number=10000):
    if p_significant_digits == 0:
        # print("this part is not finished")
        # return -1
        iterative_number = 0

        # 按照精度计算
        D = np.zeros([n, n])
        L = np.zeros([n, n])
        U = np.zeros([n, n])
        B = np.zeros([n, n])
        for i in range(n):
            for j in range(n):
                if i != j:
                    if i > j:
                        L[i, j] = -a_jacobi[i, j]
                    elif i < j:
                        U[i, j] = -a_jacobi[i, j]
                    B[i, j] = - a_jacobi[i, j] / a_jacobi[i, i]
                elif i == j:
                    D[i, j] = a_jacobi[i, j]
                    B[i, j] = 0

        value, vector = np.linalg.eig(B)
        rou_B = np.max(abs(value))
        print("||B||: " + str(rou_B))
        L = rou_B / (1 - rou_B)
        # print("L: "+str(L))
        # 利用 xk+1 冲掉 xk
        x_star_sub_xk = 0
        xk_sub_xk_1 = 0
        while x_star_sub_xk >= L * xk_sub_xk_1 and iterative_number < max_iterative_number:
            iterative_number = iterative_number + 1
            if iterative_number % 1000 == 0:
                print("迭代次数" + str(iterative_number))
            x_last_jacobi = x_jacobi.copy()
            # xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数
            for ja_i in range(n):
                ja_sum = 0
                for ja_j in range(n):
                    if ja_j != ja_i:
                        ja_sum = ja_sum + a_jacobi[ja_i, ja_j] * x_last_jacobi[ja_j, 0]

                x_jacobi[ja_i, 0] = np.longdouble((b_jacobi[ja_i, 0] - ja_sum)) / a_jacobi[ja_i, ja_i]
                # print(x_jacobi)
                # print(x_last_jacobi)
                # if iterative_number > 1:
            x_star_sub_xk = np.max(abs(x - x_jacobi))
            # print("x_star_sub_xk: " + str(x_star_sub_xk))
            xk_sub_xk_1 = np.max(abs(x_jacobi - x_last_jacobi))
            # print("xk_sub_xk_1: " + str(xk_sub_xk_1 * L))
        print("迭代次数：" + str(iterative_number))


    elif p_significant_digits == -1:
        # 直接按次数计算
        k = 0
        while k < 10000:
            k = k + 1
            print("迭代次数" + str(k))
            x_last_jacobi = x_jacobi.copy()
            # xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数
            for ja_i in range(n):
                ja_sum = 0
                for ja_j in range(n):
                    if ja_j != ja_i:
                        ja_sum = ja_sum + a_jacobi[ja_i, ja_j] * x_last_jacobi[ja_j, 0]

                x_jacobi[ja_i, 0] = (b_jacobi[ja_i, 0] - ja_sum) / a_jacobi[ja_i, ja_i]

    elif p_significant_digits == -2:
        # 按照精度计算
        k = 0
        x_star_sub_xk = 1

        while x_star_sub_xk > 0.0001:
            # x_last_jacobi = x_jacobi.copy()
            x_last_jacobi = x_jacobi.copy()
            k = k + 1
            print("迭代次数" + str(k))

            # xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数
            for ja_i in range(3):
                ja_sum = 0
                for ja_j in range(3):
                    if ja_j != ja_i:
                        ja_sum = ja_sum + a_jacobi[ja_i, ja_j] * x_last_jacobi[ja_j, 0]

                x_jacobi[ja_i, 0] = (b_jacobi[ja_i, 0] - ja_sum) / a_jacobi[ja_i, ja_i]
            x_star_sub_xk = np.max(abs(x_last_jacobi - x_jacobi))
            print("误差：" + str(x_star_sub_xk))


# jacobi_a = A.copy()
# jacobi_b = b.copy()
# jacobi_x = np.zeros([n, 1])
# jacobi_iterative_method(jacobi_a, jacobi_b, jacobi_x, 0, max_iterative_number=10000)
# np.savetxt('delta_x.txt', x - jacobi_x, fmt='%0.90f')
# print(x - jacobi_x)


"""
Gauss-Seidel
设Ax=b,其中A为n阶非奇异矩阵且 aii！= 0（i=1,2,3,...,n）
数组x(n)开始存放x0，后存放xk，N0为最大迭代次数
"""


def Gauss_Seidel_iterative_method(a_gs, b_gs, x_gs, p_significant_digits=-1, max_iterative_number=10000):
    if p_significant_digits != -1:
        # print("this part is not finished")
        # return -1
        iterative_number = 0

        # 按照精度计算
        D = np.zeros([n, n])
        L = np.zeros([n, n])
        U = np.zeros([n, n])
        B = np.zeros([n, n])
        for i in range(n):
            for j in range(n):
                if i != j:
                    if i > j:
                        L[i, j] = -a_gs[i, j]
                    elif i < j:
                        U[i, j] = -a_gs[i, j]
                    B[i, j] = - a_gs[i, j] / a_gs[i, i]
                elif i == j:
                    D[i, j] = a_gs[i, j]
                    B[i, j] = 0

        value, vector = np.linalg.eig(B)
        rou_B = np.max(abs(value))
        L = rou_B / (1 - rou_B)
        # 利用 xk+1 冲掉 xk
        while x_star_sub_xk < L * xk_sub_xk_1 or iterative_number > max_iterative_number:
            iterative_number = iterative_number + 1
            x_last_jacobi = x_gs.copy()
            # xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数
            for ja_i in range(n):
                ja_sum = 0
                for ja_j in range(n):
                    if ja_j != ja_i:
                        ja_sum = ja_sum + a_gs[ja_i, ja_j] * x_gs[ja_j, 0]

                x_gs[ja_i, 0] = np.longdouble(b_gs[ja_i, 0] - ja_sum) / a_gs[ja_i, ja_i]
                x_star_sub_xk = np.max(abs(x - x_gs))
                xk_sub_xk_1 = np.max(abs(x_gs - x_last_jacobi))
        print(iterative_number)


    elif p_significant_digits == -1:
        # 直接按次数计算
        k = 0
        while k < max_iterative_number:
            k = k + 1
            if k % 1000 == 0:
                print("迭代次数" + str(k))
            # xk+1_i = (bi- sum(j=1,j!=i)(n) aij*xk_j)/aii  i=1,2,....n  k=0,1,...表示迭代次数
            for ja_i in range(n):
                ja_sum = 0
                for ja_j in range(n):
                    if ja_j != ja_i:
                        ja_sum = ja_sum + a_gs[ja_i, ja_j] * x_gs[ja_j, 0]
                x_gs[ja_i, 0] = (b_gs[ja_i, 0] - ja_sum) / a_gs[ja_i, ja_i]

#
jacobi_a = A.copy()
jacobi_b = b.copy()
jacobi_x = np.zeros([n, 1])
Gauss_Seidel_iterative_method(jacobi_a, jacobi_b, jacobi_x, -1, max_iterative_number=10000)
np.savetxt('delta_x.txt', x - jacobi_x, fmt='%0.90f')
print(x - jacobi_x)
