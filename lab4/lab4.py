import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

C = [[8, -11, -15, -1, 0, 0],
     [0.5, 0, 0, 0, -1, 0],
     [10, 11, 15, 0, 0, -1],
     [-8, 11, 15, -1, 0, 0],
     [-0.5, 0, 0, 0, -1, 0],
     [-10, -11, -15, 0, 0, -1]]

c = [0, 0, 0, 1, 1, 1]

mid_b = [6, -9, 8, -6, 9, -8]

bounds = [[None, None] for i in range(len(mid_b))]
res_interior = opt.linprog(c, A_ub=C, b_ub=mid_b, bounds=bounds, method='interior-point')
print()
print('Interior-method results:')
print('x: (', np.around(res_interior.x[0], decimals=4), ',', np.around(res_interior.x[1], decimals=4), ',',
      np.around(res_interior.x[2], decimals=4), ')')
print('w: (', np.around(res_interior.x[3], decimals=4), ',', np.around(res_interior.x[4], decimals=4), ',',
      np.around(res_interior.x[5], decimals=4), ')')

bounds = [[None, None] for i in range(len(mid_b))]
res_simplex = opt.linprog(c, A_ub=C, b_ub=mid_b, bounds=bounds, method='simplex')
print()
print('Simplex results:')
print('x: (', np.around(res_simplex.x[0], decimals=4), ',', np.around(res_simplex.x[1], decimals=4), ',',
      np.around(res_simplex.x[2], decimals=4), ')')
print('w: (', np.around(res_simplex.x[3], decimals=4), ',', np.around(res_simplex.x[4], decimals=4), ',',
      np.around(res_simplex.x[5], decimals=4), ')')

print()
print('Delta solution for x[0]: ', abs(res_interior.x[0] - res_simplex.x[0]))



bounds_none = [[None, None] for i in range(len(mid_b))]
steps = np.arange(0, 2, 0.1)
points_2, points_3 = [], []
simplex_2, simplex_3 = [], []
for i in steps:
    bounds = ((None, None), (i, None), (i, None), (None, None), (None, None), (None, None))
    res_p = opt.linprog(c, A_ub=C, b_ub=mid_b, bounds=bounds, method='interior-point')
    res_s = opt.linprog(c, A_ub=C, b_ub=mid_b, bounds=bounds, method='simplex')
    points_2.append(res_p.x[1])
    simplex_2.append(res_s.x[1])
    points_3.append(res_p.x[2])
    simplex_3.append(res_s.x[2])
plt.figure()

plt.subplot(2, 1, 1)
plt.plot(steps, points_2, label='Interior-point')
plt.plot(steps, simplex_2, label='Simplex')
plt.legend()
plt.title('Changing of lowest bound of both 2-nd coordinate')

plt.subplot(2, 1, 2)
plt.plot(steps, points_3, label='Interior-point')
plt.plot(steps, simplex_3, label='Simplex')
plt.legend()
plt.title('Changing of lowest bound of both 3-rd coordinate')
plt.savefig('Change.png', format='png')
plt.show()

