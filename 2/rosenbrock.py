from config import *
from math import exp
import numpy as np


def test_func(x: np.array) -> float:
    """
    Test function, with params according to variant

    :param x: x-coord
    :param y: y-coord

    :return: function value
    """
    return A1 * exp(- ((x[0, 0] - a1)/(b1))**2 - ((x[0, 1] - c1)/(d1))**2) + A2 * exp(- ((x[0, 0] - a2)/(b2))**2 - ((x[0, 1] - c2)/(d2))**2)


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


def main():
    x = np.array([0, 0], dtype = float)
    s = np.matrix([[1, 0], [0, 1]], dtype = float)
    a = np.matrix([[0, 0], [0, 0]], dtype = float)
    eps = 1e-7

    for i in range(200):

        begin, end = min_interval(0, rosen_func, eps, x, s[0])
        lambda1: float = golden_section(begin, end, rosen_func, eps, x, s[0])

        x = x + s[0] * lambda1

        begin, end = min_interval(0, rosen_func, eps, x, s[1])
        lambda2: float = golden_section(begin, end, rosen_func, eps, x, s[1])

        x = x + s[1] * lambda2

        a[0] = lambda1 * s[0] + lambda2 * s[1]
        if (abs(lambda1) < abs(lambda2)):
            a[1] = lambda1 * s[1]
        else:
            a[1] = lambda2 * s[1]

        s[0] = a[0] / np.linalg.norm(a[0])

        b = a[1] - (a[1] * s[0].transpose()) * s[0]

        s[1] = b / np.linalg.norm(b)

    print(x)


if __name__ == "__main__":
    main()