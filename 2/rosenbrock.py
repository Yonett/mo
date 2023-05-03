from config import *
from math import exp, acos, sqrt
import numpy as np
import csv

funcs_count: int = 0

eps_range = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7}


def test_func(x: np.array) -> float:
    global funcs_count
    funcs_count += 1
    return A1 * exp(- ((x[0, 0] - a1)/(b1))**2 - ((x[0, 1] - c1)/(d1))**2) + A2 * exp(- ((x[0, 0] - a2)/(b2))**2 - ((x[0, 1] - c2)/(d2))**2)


def quad_func(x: np.array) -> float:
    global funcs_count
    funcs_count += 1
    return 100 * (x[0, 1] - x[0, 0])**2 + (1 - x[0, 0])**2


def rosen_func(x: np.array) -> float:
    global funcs_count
    funcs_count += 1
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


def angle(x: np.array, s: np.array):
    return acos((x[0,0]*s[0,0] + x[0,1]*s[0,1]) / (sqrt(x[0,0] ** 2 + x[0,1]**2) * sqrt(s[0,0] ** 2 + s[0,1]**2))) * 180 / 3.1415926535


def rosenbrock(x: np.array, eps: float, target: callable):
    
    s = np.matrix([[1, 0], [0, 1]], dtype = float)
    a = np.matrix([[0, 0], [0, 0]], dtype = float)

    global funcs_count
    funcs_count = 0
    iters: int = 0
    deltaF: float = eps
    deltaX: bool = True
    
    with open(f'rosenbrock/t2/{target.__name__}_{eps:1.1e}.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        while (deltaF >= eps and deltaX == True):

            iters += 1
            x_prev = np.asmatrix(x)

            begin, end = min_interval(0, target, eps, x, s[0])
            lambda1: float = golden_section(begin, end, target, eps, x, s[0])

            x = x + s[0] * lambda1

            spamwriter.writerow([iters, x[0, 0], x[0, 1], target(x), s[0, 0], s[0, 1], lambda1, abs(x[0, 0] - x_prev[0, 0]), abs(x[0, 1] - x_prev[0, 1]), abs(target(x) - target(x_prev)), angle(x, s[0])])

            begin, end = min_interval(0, target, eps, x, s[1])
            lambda2: float = golden_section(begin, end, target, eps, x, s[1])

            x = x + s[1] * lambda2

            spamwriter.writerow([iters, x[0, 0], x[0, 1], target(x), s[1, 0], s[1, 1], lambda2, abs(x[0, 0] - x_prev[0, 0]), abs(x[0,1] - x_prev[0, 1]), abs(target(x) - target(x_prev)), angle(x, s[1])])

            a[0] = lambda1 * s[0] + lambda2 * s[1]

            if (abs(lambda1) < abs(lambda2)):
                a[1] = lambda1 * s[1]
            else:
                a[1] = lambda2 * s[1]

            s[0] = a[0] / np.linalg.norm(a[0])

            b = a[1] - (a[1] * s[0].transpose()) * s[0]

            s[1] = b / np.linalg.norm(b)

            for i in range(2):
                if (abs(x[0, i] - x_prev[0, i]) < eps):
                    deltaX = False
                    break

            deltaF = abs(target(x) - target(x_prev))


        return(x, iters, funcs_count)


def main():
    x = np.array([10, -8], dtype = float)
    targets = {test_func, quad_func, rosen_func}

    for target in targets:
        with open(f'rosenbrock/t1/{target.__name__}.csv', 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for eps in eps_range:
                res, iters, funcs = rosenbrock(x, eps, target)
                spamwriter.writerow([x[0], x[1], eps, iters, funcs, res[0, 0], res[0, 1], target(res)])


if __name__ == "__main__":
    main()