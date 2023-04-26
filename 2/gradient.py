import numpy as np

def test_func(x: np.array) -> float:
    return x[0, 0]**2 + x[0, 1]**2

def quad_func(x: np.array) -> float:
    """
    Quadratic function

    :param x: vector

    :return: function value
    """
    return 100 * (x[0, 1] - x[0, 0])**2 + (1 - x[0, 0])**2


def rosen_func(x: np.array) -> float:
    """
    Rosenbrock function

    :param x: x-coord
    :param y: y-coord

    :return: function value
    """
    return 100 * (x[0, 1] - x[0, 0]**2)**2 + (1 - x[0, 0])**2


def min_interval(begin: float, func: callable, delta: float, x: np.array, s: np.array) -> tuple[float, float]:
    """
    Finds an interval for search the minimum of a function

    :param begin: start point for searching
    :param func: callable function
    :param delta: precision for first iteration

    :return: void
    """

    iters = h = 0
    prev = curr = next = begin
    if (func(x + begin * s) > func(x + (begin + delta) * s)):
        next = begin + delta
        h = delta
    if (func(x + begin * s) > func(x + (begin - delta) * s)):
        next = begin - delta
        h = -delta

    if (h != 0):
        while (func(x + curr * s) > func(x + next * s)):
            iters += 1
            prev = curr
            curr = next

            h *= 2
            next = curr + h

        return prev, next
    else:
        return begin - delta, begin + delta


def golden_section(begin: int, end: int, func: callable, eps: float, x: np.array, s: np.array) -> float:
    """
    Finds an interval that contains the minimum of a function
    Using a Golden Section Algorithm

    :param begin: begin of initial interval
    :param end: end of initial interval
    :param func: callable function (x, y)
    :param eps: precision

    :return: void
    """

    golden_number = 0.381966011

    iters = 0

    left = begin + golden_number * (end - begin)
    right = end - golden_number * (end - begin)

    value_left = func(x + left * s)
    value_right = func(x + right * s)

    while (abs(end - begin) > eps):
        iters += 1

        if (value_left < value_right):
            end = right
            right = left
            value_right = value_left
            left = begin + golden_number*(end - begin)
            value_left = func(x + left * s)
        else:
            begin = left
            left = right
            value_left = value_right
            right = end - golden_number*(end - begin)
            value_right = func(x + right * s)
    
    return (begin + end) / 2


def gradient(x: np.array, eps: float, func: callable) -> np.array:
	
    x = np.asmatrix(x)
    
    func_value = func(x)

    result = np.array([0, 0], dtype=float)
    _x = x

    for i in range(2):
        _x[0, i] += eps
        result[i] = (func(_x) - func_value) / eps
        _x[0, i] -= eps

    return result

def main():
    target: callable = rosen_func
    x = np.array([2, 2], dtype=float)
    eps: float = 1e-7

    nu = np.matrix([[1, 0], [0, 1]], dtype = float)

    grad = gradient(x, eps, target)

    for i in range(100):

        if (i % 2 == 0):
            nu = np.matrix([[1, 0], [0, 1]], dtype = float)

        if (np.linalg.det(nu) <= 0):
            nu = np.matrix([[1, 0], [0, 1]], dtype = float)

        nugrad = np.matmul(nu, grad)

        begin, end = min_interval(0, target, eps, x, nugrad)
        lambd: float = golden_section(begin, end, target, eps, x, nugrad)

        dx = lambd * nugrad
        x = x + dx

        _grad = grad
        grad = gradient(x, eps, target)

        dgrad = grad - _grad

        tmp1 = dx - np.matmul(nu, dgrad)
        tmp2 = tmp1.transpose()

        tmp3 = np.matmul(tmp2, tmp1)

        tmp4 = tmp1 * np.asmatrix(dgrad).T

        if (tmp4 == 0):
            nu = np.matrix([[1, 0], [0, 1]], dtype = float)
        else:
            dnu = tmp3 / tmp4
            nu = nu + dnu

    print(x)
    

if __name__ == "__main__":
    main()