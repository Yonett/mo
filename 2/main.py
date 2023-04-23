from config import *
from math import exp
import numpy as np


def test_func(x: int, y: int) -> float:
    """
    Test function, with params according to variant

    :param x: x-coord
    :param y: y-coord

    :return: function value
    """
    return A1 * exp(- ((x - a1)/(b1))**2 - ((y - c1)/(d1))**2) + A2 * exp(- ((x - a2)/(b2))**2 - ((y - c2)/(d2))**2)


def quad_func(x: int, y: int) -> int:
    """
    Quadratic function

    :param x: x-coord
    :param y: y-coord

    :return: function value
    """
    return 100 * (y - x)**2 + (1 - x)**2


def rosen_func(x: int, y: int) -> int:
    """
    Rosenbrock function

    :param x: x-coord
    :param y: y-coord

    :return: function value
    """
    return 100 * (y - x**2)**2 + (1 - x)**2


def rosenbrock_method(x: np.array, eps: float, func: callable):
    orthos = np.matrix([[1, 0], [0, 1]])

    

    pass


def main():
    test = np.matrix([[0, 1], [2, 3]])
    #arr = np.array([1, 2])
    arr = np.array([])
    arr = np.append(arr, 1)
    arr = np.append(arr, 2)
    print(type(np.matmul(test, arr)))
    print(np.linalg.norm(arr))
    print(test[0, 0])
    print(arr + arr * 2)


if __name__ == "__main__":
    main()