import numpy as np
import sympy as sp
from mpmath import findroot

def markBecca(midMatrix, radMatrix):
    midInvAbsMatrix = np.abs(np.linalg.inv(midMatrix))

    A = midInvAbsMatrix.dot(radMatrix)
    rho = np.max(np.abs(np.linalg.eigvals(A)))

    if rho < 1:
        print('Result Becca`s mark: Non-special matrix, rho = ', np.around(rho, decimals=3))
    else:
        print('Result Becca`s mark: Undefined, rho = ', np.around(rho, decimals=3))


def markDiagMax(midMatrix, radMatrix, size):
    midInvAbsMatrix = np.abs(np.linalg.inv(midMatrix))
    A = radMatrix.dot(midInvAbsMatrix)

    dMax = max([A[i][i] for i in range(size)])
    if dMax >= 1:
        print('Result Diagonal max mark: Special matrix, Max value in diag = ', np.around(dMax, decimals=3))
    else:
        print('Result Diagonal max mark: Undefined, Max value in diag = ', np.around(dMax, decimals=3))


def intMul(int1, int2):
    return min(int1[0] * int2[0], int1[0] * int2[1], int1[1] * int2[0], int1[1] * int2[1]), \
           max(int1[0] * int2[0], int1[0] * int2[1], int1[1] * int2[0], int1[1] * int2[1])


def intDif(int1, int2):
    return int1[0] - int2[1], int1[1] - int2[0]


def det2(matrix):
    return intDif(intMul(matrix[0][0], matrix[1][1]), intMul(matrix[0][1], matrix[1][0]))


def task1():
    eps = float(input('Task 1\nEnter eps: '))
    midMatrix = np.array([[1, 1], [1.1, 1]])
    radMatrix = np.array([[eps, eps], [eps, eps]])
    matrix = [[[1 - eps, 1 + eps], [1 - eps, 1 + eps]], [[1.1 - eps, 1.1 + eps], [1 - eps, 1 + eps]]]
    markBecca(midMatrix, radMatrix)
    markDiagMax(midMatrix, radMatrix, len(matrix))
    D = det2(matrix)
    print('det = (', np.around(D[0], decimals=3), ',', np.around(D[1], decimals=5), ')')


def task2():
    n = int(input('Task 2\nEnter dim: '))
    e = sp.Symbol("e", real=True)
    row = [e for i in range(n)]
    m = sp.Matrix([row for i in range(n)])
    for i in range(n):
        m[i, i] = 1
    d = m.det()
    print(d)

    left = [0]
    right = [0]

    for arg in d.args:
        if e not in arg.free_symbols:
            left[0] += arg
            right[0] += arg
        else:
            if arg.is_nonpositive:
                left[0] += arg
            else:
                right[0] += arg

    ans = findroot(lambda x: left[0].evalf(subs={e: x}), 0)
    print('Left bound: ', left)
    print('Right bound: ', right)
    print('Determinant: (', left[0].subs(e, ans), ', ', right[0].subs(e, ans), ')')
    print('Epsilon:', ans)
    print('1/(N-1) = ', 1. / (n - 1))


#task1()
task2()