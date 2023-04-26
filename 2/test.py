import numpy as np

x = np.array([1, 1], dtype=float)
nu = np.matrix([[1, 0], [0, 1]], dtype = float)
grad = np.array([2, 2], dtype=float)
_grad = np.array([-6, -6], dtype=float)

lambd = 2

nugrad = np.matmul(nu, grad)

print(nugrad)

dx = lambd * nugrad

print(dx)

x = x - dx

print(x)

dx = -dx

dg = _grad - grad

print(dg)

tmp1 = dx - np.matmul(nu, dg)
tmp2 = tmp1.transpose()

print(tmp1)
print(tmp2)
tmp3 = np.matmul(tmp2, tmp1)
print(tmp3)

tmp4 = tmp1 * np.asmatrix(dg).T

dnu = tmp3 / tmp4

print(dnu)

# _grad = grad
# grad = gradient(x, eps, target)

# dgrad = grad - _grad

# temp = dx - np.matmul(nu, grad)
# temp = np.ravel(temp)

# dnu = np.matmul(temp, temp) / np.matmul(temp, dgrad)

# if (np.matmul(temp, dgrad) == 0):
#     nu = np.matrix([[1, 0], [0, 1]], dtype = float)
# else:
#     nu = nu + dnu

# print(x)