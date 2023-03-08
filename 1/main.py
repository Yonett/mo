import matplotlib.pyplot as plt
import numpy as np
import csv

eps_range = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
dichotomy_objectives_range = []
golden_section_objectives_range = []
fibonacci_objectives_range = []

def objective(func, x):
    return eval(func, {"x": x})


def fibonacci_number(n: int):
    prev = curr = 1
    n -= 2

    while (n > 0):
        prev, curr = curr, prev + curr
        n -= 1

    return curr


def dichotomy(begin: float, end: float, func: str, eps: float, delta: float):
    """
    Finds an interval that contains the minimum of a function
    Using Dichotomy Algorithm
    Also creates a file with csv table of the following format:

    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+
    | iter | left | right | value_left | value_right | begin(i) | end(i) | end(i) - begin(i) | end(i - 1) - begin(i - 1) /  end(i) - begin(i) |
    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+

    File with table has a name "dichotomy_{eps:1.1e}.csv"

    :param begin: begin of initial interval
    :param end: end of initial interval
    :param func: function expression in Python format, must be
                 valid for standard eval() function
    :param eps: precision
    :param delta: value of deflection from the middle of interval

    :return: void
    """

    iters = 0
    prev_end_begin = 0
    func_string = func.replace("**", "^").replace("\n", "")

    with open(f'dichotomy\{func_string}_{eps:1.1e}.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        while (abs(end - begin) > eps):

            iters += 1

            left = (begin + end - delta) / 2
            right = (begin + end + delta) / 2

            value_left = objective(func, left)
            value_right = objective(func, right)

            spamwriter.writerow([iters, left, right, value_left, value_right, begin, end, end - begin, prev_end_begin / (end - begin)])

            prev_end_begin = end - begin

            if (value_left < value_right):
                end = right
            elif (value_right < value_left):
                begin = left
            else:
                begin = left
                end = right

    dichotomy_objectives_range.append(iters * 2)
    print(f"{begin:.8e}, {end:.8e}, {iters}")


def golden_section(begin: float, end: float, func: str, eps: float):
    """
    Finds an interval that contains the minimum of a function
    Using a Golden Section Algorithm
    Also creates a file with csv table of the following format:

    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+
    | iter | left | right | value_left | value_right | begin(i) | end(i) | end(i) - begin(i) | end(i - 1) - begin(i - 1) /  end(i) - begin(i) |
    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+

    File with table has a name "golden_section_{eps:1.1e}.csv"

    :param begin: begin of initial interval
    :param end: end of initial interval
    :param func: function expression in Python format, must be
                 valid for standard eval() function
    :param eps: precision

    :return: void
    """

    golden_number = 0.381966011

    iters = 0
    prev_end_begin = 0
    func_string = func.replace("**", "^").replace("\n", "")

    left = begin + golden_number * (end - begin)
    right = end - golden_number * (end - begin)

    value_left = objective(func, left)
    value_right = objective(func, right)

    with open(f'golden_section\{func_string}_{eps:1.1e}.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        while (abs(end - begin) > eps):
            iters += 1

            spamwriter.writerow([iters, left, right, value_left, value_right, begin, end, end - begin, prev_end_begin / (end - begin)])

            prev_end_begin = end - begin

            if (value_left < value_right):
                end = right
                right = left
                value_right = value_left
                left = begin + golden_number*(end - begin)
                value_left = objective(func, left)
            else:
                begin = left
                left = right
                value_left = value_right
                right = end - golden_number*(end - begin)
                value_right = objective(func, right)

    golden_section_objectives_range.append(iters + 2)
    print(f"{begin:.8e}, {end:.8e}, {iters}")


def fibonacci(begin: float, end: float, func: str, eps: float):
    """
    Finds an interval that contains the minimum of a function
    Using a Fibonacci Algorithm
    Also creates a file with csv table of the following format:

    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+
    | iter | left | right | value_left | value_right | begin(i) | end(i) | end(i) - begin(i) | end(i - 1) - begin(i - 1) /  end(i) - begin(i) |
    +------+------+-------+------------+-------------+----------+--------+-------------------+------------------------------------------------+

    File with table has a name "fibonacci_{eps:1.1e}.csv"

    :param begin: begin of initial interval
    :param end: end of initial interval
    :param func: function expression in Python format, must be
                 valid for standard eval() function
    :param eps: precision

    :return: void
    """
    n = 0
    iters = 0
    prev_end_begin = 0
    func_string = func.replace("**", "^").replace("\n", "")

    while (fibonacci_number(n + 2) < ((end - begin) / eps)):
        n +=1

    left = begin + (fibonacci_number(n) * (end - begin) / fibonacci_number(n + 2))
    right = begin + (fibonacci_number(n + 1) * (end - begin) / fibonacci_number(n + 2))

    value_left = objective(func, left)
    value_right = objective(func, right)

    with open(f'fibonacci\{func_string}_{eps:1.1e}.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        while (abs(end - begin) > eps):
            iters += 1

            spamwriter.writerow([iters, left, right, value_left, value_right, begin, end, end - begin, prev_end_begin / (end - begin)])

            prev_end_begin = end - begin

            if (value_left < value_right):
                end = right

                right = left

                value_right = value_left
                left = begin + (fibonacci_number(n - iters + 1) * (end - begin) / fibonacci_number(n - iters + 3))
                value_left = objective(func, left)

            else:
                begin = left

                left = right
                value_left = value_right

                right = begin + (fibonacci_number(n - iters + 2) * (end - begin) / fibonacci_number(n - iters + 3))
                value_right = objective(func, right)


    fibonacci_objectives_range.append(iters + 2)
    print(f"{begin:.8e}, {end:.8e}, {iters}")


def min_interval(begin: float, func: str, delta: float):
    iters = h = 0
    prev = curr = next = begin

    if (objective(func, begin) > objective(func, begin + delta)):
        next = begin + delta
        h = delta
    if (objective(func, begin) > objective(func, begin - delta)):
        next = begin - delta
        h = -delta

    if (h != 0):
        while (objective(func, curr) > objective(func, next)):
            iters += 1
            prev = curr
            curr = next

            h *= 2
            next = curr + h

        print(f"[{prev},{next}], iters = {iters}")
    else:
        print(f"[{begin - delta}, {begin + delta}], iters = 0")


def main():
    file = open("funcs.txt")
    for func in file:

        dichotomy_objectives_range.clear()
        golden_section_objectives_range.clear()
        fibonacci_objectives_range.clear()

        for eps in eps_range:

            delta = eps * 1.e-3
            print(f"{func}\nprecision: {eps}\tdelta: {delta}")

            print("Dichotomy:")
            dichotomy(-2, 20, func, eps, delta)

            print("Golden Section:")
            golden_section(-2, 20, func, eps)

            print("Fibonacci:")
            fibonacci(-2, 20, func, eps)

        print(f"{func}delta: {delta}")
        print("Min interval (begin = 0.0):")
        min_interval(0.0, func, delta)
        print("-------------------------")

    fig, ax = plt.subplots()
    ax.plot(np.log10(eps_range), dichotomy_objectives_range)
    ax.plot(np.log10(eps_range), golden_section_objectives_range)
    ax.plot(np.log10(eps_range), fibonacci_objectives_range)
    ax.grid()
    plt.show()
    
    file.close()


if __name__ == "__main__":
    main()