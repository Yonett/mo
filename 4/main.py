from config import *
from math import log, ceil
import random
import numpy as np
import csv

eps_range = [1, 1e-1, 1e-2]
P_range = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
func_count = 0

def func(x: np.array):
    global func_count
    func_count += 1
    result = 0
    for i in range(6):
        result += C[i] / (1+(x[0, 0] - a[i])**2 + (x[0, 1] - b[i])**2)
    return -result


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


def simple_random_search(P: float, eps: float, target: callable):
    dot = np.array([], dtype=float)
    V: float = 400.0

    V_eps: float = eps**2
    P_eps = V_eps / V

    N: float = log((1 - P))/log(1 - P_eps)
    N: int = ceil(N)

    min = np.array([], dtype=float)
    value = 999_999

    i = 0
    while i < (N * 2):
        random.seed(i)
        dot = np.append(dot, random.uniform(-10, 10))

        if i % 2 != 0:
            if (target(np.asmatrix(dot)) < value):
                value = target(np.asmatrix(dot))
                min = dot
            dot = np.array([], dtype=float)

        i += 1

    return min, -value, N

def make_srs_research() -> None:
    with open(f'simple_random_search.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for eps in eps_range:
            for P in P_range:
                dot, value, N = simple_random_search(P, eps, func)
                spamwriter.writerow([eps, P, N, dot[0], dot[1], value])

def alg1(seed: float, m: int, eps: float, target: callable):
    random.seed(seed)
    dot = np.array([], dtype=float)
    dot = np.append(dot, random.uniform(-10, 10))
    dot = np.append(dot, random.uniform(-10, 10))

    result = rosenbrock(dot, eps, target)
    _result = result

    attempts = 0

    while attempts < m:
        seed += 1
        random.seed(seed)
        dot = np.array([], dtype=float)
        dot = np.append(dot, random.uniform(-10, 10))
        dot = np.append(dot, random.uniform(-10, 10))

        result = rosenbrock(dot, eps, target)

        if (target(result) < target(_result)):
            attempts = 0
            _result = result
        else:
            attempts += 1

    return _result, -target(_result)

def alg2(seed: float, m: int, eps: float, target: callable):
    random.seed(seed)
    dot = np.array([], dtype=float)
    dot = np.append(dot, random.uniform(-10, 10))
    dot = np.append(dot, random.uniform(-10, 10))

    result = rosenbrock(dot, eps, target)
    value = target(result)

    attempts = 0

    while (attempts < m):
        seed += 1
        random.seed(seed)
        dot = np.array([], dtype=float)
        dot = np.append(dot, random.uniform(-10, 10))
        dot = np.append(dot, random.uniform(-10, 10))

        if (target(np.asmatrix(dot)) < value):
            attempts = 0
            result = rosenbrock(dot, eps, target)
            value = target(result)
        else:
            attempts += 1

    return result, -value
   
def alg3(seed: float, m: int, eps: float, target: callable):
    random.seed(seed)
    dot = np.array([], dtype=float)
    dot = np.append(dot, random.uniform(-10, 10))
    dot = np.append(dot, random.uniform(-10, 10))

    result = rosenbrock(dot, eps, target)
    value = target(result)

    attempts = 0

    while (attempts < m):
        seed += 1
        rand_dir = np.array([], dtype=float)
        rand_dir = np.append(rand_dir, random.uniform(-10, 10))
        rand_dir = np.append(rand_dir, random.uniform(-10, 10))

        step: int = 1
        temp_dot = dot + (dot - rand_dir) * step
        temp_func = target(np.asmatrix(temp_dot))

        while (temp_func >= value and
               temp_dot[0] > -10 and temp_dot[0] < 10 and
               temp_dot[1] > -10 and temp_dot[1] < 10):
            
            attempts += 1
            step += 1
            temp_dot = dot + (dot - rand_dir) * step
            temp_func = target(np.asmatrix(temp_dot))
        
        if (temp_func < value):
            value = temp_func
            result = temp_dot
            attempts = 0


    return result, -value

def main():
    #make_srs_research()
    algs = [alg1, alg2, alg3]
    ms = [10, 100, 1000, 10000]
    seeds = [1, 10, 100, 1000, 10000]
    global func_count

    for seed in seeds:
        with open(f'seed_{seed}.csv', 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for m in ms:
                for alg in algs:
                    print(f"{seed}, {m}, {alg.__name__} working...")
                    func_count = 0
                    dot, value = alg(seed, m, 1e-7, func)
                    dot = np.asmatrix(dot)
                    spamwriter.writerow([m, alg.__name__, func_count, dot[0, 0], dot[0, 1], value])
                    print(f"{seed}, {m}, {alg.__name__} done")

    # print(alg1(1, 10, 1e-7, func))
    # print(alg2(1, 10, 1e-7, func))
    # print(alg3(1, 10, 1e-7, func))

if __name__ == "__main__":
    main()