import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random


def createHeatMap(matrix, title='Matrix'):
    plt.figure(figsize=(16, 9))
    sns.heatmap(matrix, center=0)
    plt.title(title)
    plt.savefig(title + '.png', format='png')
    plt.show()


def plot_graph(x, x_inf, x_sup, name='Picture.png'):
    plt.figure(figsize=(16, 9))
    plt.plot(x, label='Исходное сгенерированное решение')
    plt.plot(x_inf, label='Inf')
    plt.plot(x_sup, label='Sup')

    plt.xlabel('X[i]')
    plt.ylabel('Value')
    plt.title('Сравнение решения с исходными данными')
    plt.legend()
    plt.savefig(name + '.png', format='png')
    plt.show()


def plot_graph2(x, f_x_inf, f_x_sup, s_x_inf, s_x_sup, cross_x_inf, cross_x_sup, name='Picture'):
    plt.figure(figsize=(16, 9))
    plt.plot(x, label='Исходное сгенерированное решение')

    plt.plot(f_x_inf, label='Inf для решения 1 матрицы')
    plt.plot(f_x_sup, label='Sup для решения 1 матрицы')

    plt.plot(s_x_inf, 'b', label='Inf для решения 2 матрицы', ls='--')
    plt.plot(s_x_sup, 'r', label='Sup для решения 2 матрицы', ls='--')

    plt.plot(cross_x_inf, 'bo', label='Inf для пересечения решений')
    plt.plot(cross_x_sup, 'ro', label='Sup для пересечения решений')

    plt.xlabel('X[i]')
    plt.ylabel('Value')
    plt.title('Сравнение решения с исходными данными')
    plt.legend()
    plt.savefig(name + '.png', format='png')
    plt.show()


def findBestMatrix(matrix, size):
    print('Task ', size, 'x', size)
    det = 0
    out_matr = []
    indexes = [i for i in range(matrix.shape[0])]
    count_iter = 500000
    for i in range(count_iter):
        cur_indexes = random.sample(indexes, size)
        cur_matr = []
        for elem in cur_indexes:
            if len(cur_matr) == 0:
                cur_matr = np.array([matrix[elem]])
            else:
                cur_matr = np.append(cur_matr, [matrix[elem]], axis=0)

        cur_det = np.linalg.det(cur_matr)
        if cur_det > det:
            det = cur_det
            out_matr = cur_matr.copy()
        if i % (count_iter / 10) == 0:
            print('Brute force progress: ', i / (count_iter / 100), '%')
            print('Current max det: ', det)

    print()
    print('Max det in find matrix :', det)
    print()
    if len(out_matr) != 0:
        return out_matr
    else:
        return matrix[:size][:size]


def createSignBlockMatrix(matrix):
    pos, neg = matrix.copy(), matrix.copy()
    pos[pos < 0] = 0
    neg[neg > 0] = 0
    neg = np.abs(neg)
    return np.block([[pos, neg], [neg, pos]])


def solve(A, b_inf, b_sup):
    sti_vec = np.append(-b_inf, b_sup)
    A_block = createSignBlockMatrix(A)

    res = np.dot(np.linalg.inv(A_block), sti_vec)

    middle_index = int(res.shape[0] / 2)
    x_inf, x_sup = -res[:middle_index], res[middle_index:]

    return x_inf, x_sup


def makeAbsoluteNonSpecial(matrix):
    if abs(np.linalg.det(matrix)) >= 1e-2 and abs(np.linalg.det(np.abs(matrix))) >= 1e-2:
        return matrix, False
    else:
        D = np.diag(np.abs(matrix))
        E = np.sum(np.abs(matrix), axis=1) - D
        for i in range(matrix.shape[0]):
            if D[i] - E[i] <= 0:
                matrix[i][i] += E[i] - D[i]
        return matrix, True


def generateRightPart(matrix):
    n = matrix.shape[0]
    x = np.random.uniform(size=n)
    rad = np.random.uniform(low=0.1, high=1.5, size=n)

    b = np.dot(matrix, x)
    b_inf = b - rad
    b_sup = b + rad

    return b_inf, b_sup, x


def calcD(D, i, j, cur_interval, b_inf, b_sup):
    size = int(D.shape[0] / 2)

    a_inf = cur_interval[0]
    a_sup = cur_interval[1]

    k = 0
    m = 0
    if a_inf * a_sup > 0:
        k = 0 if a_inf > 0 else 2
    else:
        k = 1 if a_inf < a_sup else 3

    if b_inf * b_sup > 0:
        m = 1 if b_inf > 0 else 3
    else:
        m = 2 if b_inf <= b_sup else 4

    case = 4 * k + m
    if case == 1:
        D[i][j] = a_inf
        D[i + size][j + size] = a_sup
    elif case == 2:
        D[i][j] = a_sup
        D[i + size][j + size] = a_sup
    elif case == 3:
        D[i][j] = a_sup
        D[i + size][j + size] = a_inf
    elif case == 4:
        D[i][j] = a_inf
        D[i + size][j + size] = a_inf
    elif case == 5:
        D[i][j + size] = a_inf
        D[i + size][j + size] = a_sup
    elif case == 6:
        if a_inf * b_sup < a_sup * b_inf:
            D[i][j + size] = a_inf
        else:
            D[i][j] = a_sup
        if a_inf * b_inf > a_sup * b_sup:
            D[i + size][j] = a_inf
        else:
            D[i + size][j + size] = a_sup
    elif case == 7:
        D[i][j] = a_sup
        D[i + size][j] = a_inf
    elif case == 8:
        return D
    elif case == 9:
        D[i][j + size] = a_inf
        D[i + size][j] = a_sup
    elif case == 10:
        D[i][j + size] = a_inf
        D[i + size][j] = a_inf
    elif case == 11:
        D[i][j + size] = a_sup
        D[i + size][j] = a_inf
    elif case == 12:
        D[i][j + size] = a_sup
        D[i + size][j] = a_sup
    elif case == 13:
        D[i][j] = a_inf
        D[i + size][j] = a_sup
    elif case == 14:
        return D
    elif case == 15:
        D[i][j + size] = a_sup
        D[i + size][j + size] = a_inf
    elif case == 16:
        if a_inf * b_inf > a_sup * b_sup:
            D[i][j] = a_inf
        else:
            D[i][j + size] = -a_sup
        if a_inf * b_sup < a_sup * b_inf:
            D[i + size][j + size] = a_inf
        else:
            D[i + size][j] = a_sup
    return D


def calcG(A, last_x, sti_vec):
    half_size = int(last_x.shape[0] / 2)
    inf, sup = -last_x[:half_size], last_x[half_size:]
    intervals = np.array([[inf[i], sup[i]] for i in range(inf.shape[0])])
    dot = [sum([A[i][j] * intervals[j] for j in range(len(intervals))]) for i in range(len(A))]
    dot_lower, dot_upper = np.array([i[0] for i in dot]), np.array([i[1] for i in dot])
    return np.append(-dot_lower, dot_upper) - sti_vec


def crossSolution(first_sol_x, second_sol_x):
    cross = []
    middle_index = int(first_sol_x.shape[0] / 2)
    f_inf, f_sup = -first_sol_x[:middle_index], first_sol_x[middle_index:]
    first_sol_x = np.array([f_inf, f_sup]).T
    c_inf, c_sup = -second_sol_x[:middle_index], second_sol_x[middle_index:]
    second_sol_x = np.array([c_inf, c_sup]).T

    dual_f = [False for i in range(len(first_sol_x))]
    dual_s= [False for i in range(len(first_sol_x))]

    for i in range(len(first_sol_x)):
        if first_sol_x[i][0] >= first_sol_x[i][1]:
            dual_f[i] = True
            first_sol_x[i].sort()

        if second_sol_x[i][0] >= second_sol_x[i][1]:
            dual_s[i] = True
            second_sol_x[i].sort()

        if second_sol_x[i][1] < first_sol_x[i][0] or first_sol_x[i][1] < second_sol_x[i][0]:
            cross.append([-1e+200, 1e+200])
        elif second_sol_x[i][0] <= first_sol_x[i][0] <= second_sol_x[i][1] <= first_sol_x[i][1]:
            cross.append([first_sol_x[i][0], second_sol_x[i][1]])
        elif second_sol_x[i][0] <= first_sol_x[i][0] and first_sol_x[i][1] <= second_sol_x[i][1]:
            cross.append(first_sol_x[i])
        elif first_sol_x[i][0] <= second_sol_x[i][0] <= first_sol_x[i][1] <= second_sol_x[i][1]:
            cross.append([second_sol_x[i][0], first_sol_x[i][1]])
        elif first_sol_x[i][0] <= second_sol_x[i][0] and second_sol_x[i][1] <= first_sol_x[i][1]:
            cross.append(second_sol_x[i])

        if dual_f[i] or dual_s[i]:
            tmp = cross[i][0]
            cross[i][0] = cross[i][1]
            cross[i][1] = tmp

    cross = np.array(cross)
    return f_inf, f_sup, c_inf, c_sup, cross.T[:][0], cross.T[:][1]


def subdiff2(A_inf, A_sup, b_inf, b_sup, lr=0.5, acc=1e-10, max_iter=50):
    size = A_inf.shape[0]

    A = np.array([[[A_inf[i][j], A_sup[i][j]] for j in range(size)] for i in range(size)])
    A_mid = np.array([[(A_inf[i][j] + A_sup[i][j]) / 2 for j in range(size)] for i in range(size)])

    A_block = createSignBlockMatrix(A_mid)
    sti_vec = np.append(-b_inf, b_sup)
    cur_x = np.dot(np.linalg.inv(A_block), sti_vec)

    last_x = cur_x

    iter = 0
    firstStep = True
    while firstStep or np.linalg.norm(cur_x - last_x) > acc:
        firstStep = False
        iter += 1
        if iter > max_iter:
            print('Too many iterations')
            break
        last_x = cur_x
        D = np.zeros([2 * size, 2 * size])
        for i in range(size):
            for j in range(size):
                cur_interval = A[i][j]
                last_inf = -last_x[j]
                last_sup = last_x[j + size]
                D = calcD(D, i, j, cur_interval, last_inf, last_sup)
        G = calcG(A, last_x, sti_vec)
        dx = np.dot(np.linalg.inv(D), -G)
        cur_x = last_x + lr * dx
    return cur_x, iter


def task1():
    print('------------------------------Task 1 start------------------------------------')
    A_inf = np.array([[1, 1],
                      [0, 0.1]])

    A_sup = A_inf.copy()

    b_inf, b_sup = np.array([2.1, 0.5]), np.array([1.9, 0.8])

    sol_x, iter = subdiff2(A_inf, A_sup, b_inf, b_sup)

    middle_index = int(sol_x.shape[0] / 2)
    c_inf, c_sup = -sol_x[:middle_index], sol_x[middle_index:]

    print('Iterations: ', iter)
    print('Solution: \n', np.block([[c_inf], [c_sup]]).T)
    print('------------------------------Task 1 finish------------------------------------')


def task2():
    print('------------------------------Task 2 start------------------------------------')
    file = 'matrix_n_phi_1.txt'
    A = np.loadtxt(file)

    size = min(A.shape)
    x = np.random.uniform(size=size)
    rad = np.random.uniform(low=0.1, high=1.5, size=size)

    print('---------------First matrix solution----------------')

    A_inf = findBestMatrix(A, size)
    A_sup = A_inf.copy()

    createHeatMap(A_inf, title=file + ', ' + str(A_inf.shape[0]) + 'x' + str(A_inf.shape[0]) + '_1')

    b = np.dot(A_inf, x)
    b_inf = b - rad
    b_sup = b + rad

    first_sol_x, iter = subdiff2(A_inf, A_sup, b_inf, b_sup)

    print('Iterations count for first matrix: ', iter)

    print('---------------Second matrix solution----------------')
    A_inf = findBestMatrix(A, size)
    A_sup = A_inf.copy()

    createHeatMap(A_inf, title=file + ', ' + str(A_inf.shape[0]) + 'x' + str(A_inf.shape[0]) + '_2')

    b = np.dot(A_inf, x)
    b_inf = b - rad
    b_sup = b + rad

    second_sol_x, iter = subdiff2(A_inf, A_sup, b_inf, b_sup)

    print('Iterations count for second matrix: ', iter)

    f_x_inf, f_x_sup, s_x_inf, s_x_sup, cross_x_inf, cross_x_sup = crossSolution(first_sol_x, second_sol_x)

    plot_graph2(x, f_x_inf, f_x_sup, s_x_inf, s_x_sup, cross_x_inf, cross_x_sup, name='GraphForTask2')
    print('------------------------------Task 2 finish------------------------------------')


def task3():
    print('------------------------------Task 3 start------------------------------------')
    file = 'matrix_n_phi_6.txt'
    A = np.loadtxt(file)
    A = A[:min(A.shape)][:min(A.shape)]

    createHeatMap(A, title=file + ', ' + str(A.shape[0]) + 'x' + str(A.shape[0]))
    A, needed_change = makeAbsoluteNonSpecial(A)
    if needed_change:
        createHeatMap(A, title=file + ', ' + str(A.shape[0]) + 'x' + str(A.shape[0]) + ', Revised')

    b_inf, b_sup, x = generateRightPart(A)

    sol_x, iter = subdiff2(A, A.copy(), b_inf, b_sup)

    middle_index = int(sol_x.shape[0] / 2)
    x_inf, x_sup = -sol_x[:middle_index], sol_x[middle_index:]

    print('Iterations count: ', iter)
    plot_graph(x, x_inf, x_sup, name='GraphForTask3')
    print('------------------------------Task 3 finish------------------------------------')


task1()
task2()
task3()
