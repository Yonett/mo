import numpy as np

s = np.matrix([[1, 0], [0, 1]], dtype = float)
a = np.matrix([[0, 0], [0, 0]], dtype = float)

lambda1 = 3.0
lambda2 = 4.0

a[0] = lambda1 * s[0] + lambda2 * s[1]

if (abs(lambda1) < abs(lambda2)):
    a[1] = lambda1 * s[1]
else:
    a[1] = lambda2 * s[1]

print(a)

s[0] = a[0] / np.linalg.norm(a[0])

print(s)

b = a[1] - (a[1] * s[0].transpose()) * s[0]

s[1] = b / np.linalg.norm(b)