import numpy as np
from math import log

alpha: int = 4
r: float = 1

def func_f(x: np.array) -> float:
    return 4*(x[0, 0] + x[0, 1])**2 + x[0, 0]**2

def func_g(x: np.array) -> float:
    return x[0, 0] - x[0, 1] + 5

def func_h(x: np.array) -> float:
    return x[0, 0] - 2 + x[0, 1]

def H1(x: np.array) -> float:
    return abs(func_h(x))

def H2(x: np.array) -> float:
    return abs(func_h(x)) ** 2

def H3(x: np.array) -> float:
    return abs(func_h(x)) ** alpha

def G1(x: np.array) -> float:
    value = func_g(x)
    return (value + abs(value)) / 2

def G2(x: np.array) -> float:
    value = func_g(x)
    return ((value + abs(value)) / 2) ** 2

def G3(x: np.array) -> float:
    value = func_g(x)
    return ((value + abs(value)) / 2) ** alpha

def G4(x: np.array) -> float:
    return -(1 / func_g(x))

def G5(x: np.array) -> float:
    return -log(-func_g(x))

#Функции для метода штрафный функций

def Q_0_1(x: np.array) -> float:
    return(func_f(x) + r * (G1(x)))

def Q_0_2(x: np.array) -> float:
    return(func_f(x) + r * (G2(x)))

def Q_0_3(x: np.array) -> float:
    return(func_f(x) + r * (G3(x)))

def Q_1_0(x: np.array) -> float:
    return(func_f(x) + r * (H1(x)))

def Q_2_0(x: np.array) -> float:
    return(func_f(x) + r * (H2(x)))

def Q_3_0(x: np.array) -> float:
    return(func_f(x) + r * (H3(x)))

##########################

#Функции для метода барьерных функций

def Q_0_4(x: np.array) -> float:
    return(func_f(x) + r * (G4(x)))

def Q_0_5(x: np.array) -> float:
    return(func_f(x) + r * (G5(x)))

##########################

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


def rosenbrock(x: np.array, eps: float, target: callable):
    
    s = np.matrix([[1, 0], [0, 1]], dtype = float)
    a = np.matrix([[0, 0], [0, 0]], dtype = float)

    deltaF: float = eps
    deltaX: bool = True
    
    while (deltaF >= eps and deltaX == True):

        x_prev = x

        begin, end = min_interval(0, target, eps, x, s[0])
        lambda1: float = golden_section(begin, end, target, eps, x, s[0])

        x = x + s[0] * lambda1

        begin, end = min_interval(0, target, eps, x, s[1])
        lambda2: float = golden_section(begin, end, target, eps, x, s[1])

        x = x + s[1] * lambda2

        a[0] = lambda1 * s[0] + lambda2 * s[1]

        if (abs(lambda1) < abs(lambda2)):
            a[1] = lambda1 * s[1]
        else:
            a[1] = lambda2 * s[1]

        s[0] = a[0] / np.linalg.norm(a[0])

        b = a[1] - (a[1] * s[0].transpose()) * s[0]

        s[1] = b / np.linalg.norm(b)

        for i in range(2):
            if (abs(np.asmatrix(x)[0, i] - np.asmatrix(x_prev)[0, i]) < eps):
                deltaX = False
                break

        deltaF = abs(target(np.asmatrix(x)) - target(np.asmatrix(x_prev)))

    return x

##########################

def main():

    # A solve:
    # 5,8823529411764705882352941176471
    # B solve:
    # 16

    A_solve = 5.8823529411764705882352941176471
    B_solve = 16.0

    x = np.array([0, 0], dtype = float)
    eps: float = 1e-7


    target: callable = Q_0_3
    global r
    iters = 0
    solve = rosenbrock(x, eps, target)

    while (func_g(solve) > 0):
        r *= 2
        iters += 1
        solve = rosenbrock(x, eps, target)
    print(solve, target(solve), r, iters, target(solve) - A_solve)

    r = 1
    x = np.array([0, 0], dtype = float)
    target: callable = Q_1_0

    solve = rosenbrock(x, eps, target)

    while (abs(func_h(solve)) > eps):
        r *= 2
        solve = rosenbrock(x, eps, target)
    print(solve, target(solve), r, target(solve) - B_solve)

if __name__ == "__main__":
    main()