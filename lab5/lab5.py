import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


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


def createSignBlockMatrix(matrix):
    pos, neg = matrix.copy(), matrix.copy()
    pos[pos < 0] = 0
    neg[neg > 0] = 0
    neg = np.abs(neg)
    return np.block([[pos, neg], [neg, pos]])


def solve(A, b_inf, b_sup):
    sti_vec = np.append(-b_inf, b_sup)
    A_block = createSignBlockMatrix(A)
    print(np.linalg.inv(A_block))
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


def task1():
    A = np.array([[1, 1],
                  [0, 0.1]])

    b_inf, b_sup = np.array([2.1, 0.5]), np.array([1.9, 0.8])

    x_inf, x_sup = solve(A, b_inf, b_sup)

    print(np.block([[x_inf], [x_sup]]).T)


def task2():
    files = ['matrix_n_phi_1.txt', 'matrix_n_phi_6.txt']
    for file in files:
        A = np.loadtxt(file)
        A = A[:min(A.shape)][:min(A.shape)]

        createHeatMap(A, title=file + ', ' + str(A.shape[0]) + 'x' + str(A.shape[0]))
        A, needed_change = makeAbsoluteNonSpecial(A)
        if needed_change:
            createHeatMap(A, title=file + ', ' + str(A.shape[0]) + 'x' + str(A.shape[0]) + ', Revised')

        b_inf, b_sup, x = generateRightPart(A)

        x_inf, x_sup = solve(A, b_inf, b_sup)

        plot_graph(x, x_inf, x_sup, name=file[:-4])


task1()
task2()
